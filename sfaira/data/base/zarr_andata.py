import anndata
import dask.array
import dask.dataframe
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Union
import zarr


def write_zarr(store: Union[str, Path], adata, chunks, compression_kwargs):
    """
    Prepare zarr-dask optimised on disk representation of AnnData instance.

    A zarr group (written by anndata._io.write_zarr) of:
        - .X: as zarr array which can be interfaced with dask
        - .obs: as zarr array which can be interfaced with dask but is effectively not used during reading
        - .var: as zarr array which can be interfaced with dask
        - .uns: as zarr array which can be interfaced with dask
    This group can be easily interfaced via anndata.read_zarr. Additionally, also the following representations are
    saved:

    A parquet of
        - .obs which can be interfaced with pandas.DataFrame (and dask.dataframe.DataFrame)

    :param store: File name of the store (zarr group).
    :param adata: Anndata to save.
    :param chunks: Chunking of .X for zarr.
    :param compression_kwargs: Compression kwargs for zarr.
    """
    # The write_zarr method of AnnData object does not propagate kwargs yet, use raw function here:
    anndata._io.write_zarr(store=store, adata=adata, chunks=chunks, **compression_kwargs)
    # Also write .obs as a separate file as this can be easily interfaced with DataFrames.
    # The .obs written by write_zarr is a zarr group which is evaluated when the Anndata object is initialised,
    # so cannot be kept out of memory and does not allow easy column subsetting during reading.
    adata.obs.to_parquet(path=store + "_obs.parquet", engine='pyarrow', compression='snappy', index=None)


def read_zarr(store: Union[str, Path], use_dask: bool = True, columns: Union[None, List[str]] = None) -> \
        Tuple[anndata.AnnData, Union[dask.dataframe.DataFrame, pd.DataFrame]]:
    """
    Assembles a reduced AnnData instance based on a zarr store written with anndata.

    The anndata instance inclused .obs, .var, .uns and .X.
    In contrast to anndata._io.zarr.read_zarr(), this functon interfaces components of anndata as lazy dask arrays on
    the zarr store rather than loading the arrays from the zarr store fully into memory.
    In particular, if use_dask is True:

        - .X is interfaced as a dask Array
    TODO in the future potentially:
        - .obs is interfaced from a parquet file through pyarrow as a dask DataFrame

    :param store: Path to zarr group.
    :param use_dask: Whether to use lazy dask arrays where appropriate.
    :param columns: Which columns to read into the obs copy in the output, see pandas.read_parquet().
    :return: Tuple of:
        - AnnData with .X as dask array.
        - obs table separately as dataframe
    """
    from anndata._io.zarr import _clean_uns, read_dataframe, read_attribute
    if isinstance(store, Path):
        store = str(store)

    f = zarr.open(store, mode="r")
    if use_dask:
        x = dask.array.from_zarr(url=store, component="X")
    else:
        x = read_attribute(f["X"])
    # Could use dask DataFrame here too:
    obs = pd.read_parquet(store + "_obs.parquet", columns=columns, engine="pyarrow")
    # Convert to categorical variables where possible to save memory:
    for k, dtype in zip(list(obs.columns), obs.dtypes):
        if dtype == "object":
            obs[k] = obs[k].astype(dtype="category")
    d = {
        "var": read_dataframe(f["var"]),
        "uns": read_attribute(f["uns"]),
    }
    _clean_uns(d)
    # Assemble AnnData without obs to save memory:
    adata = anndata.AnnData(**d, shape=x.shape)
    # Need to add these attributes after initialisation so that they are not evaluated:
    adata.X = x
    return adata, obs
