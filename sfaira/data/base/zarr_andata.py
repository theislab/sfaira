import anndata
import dask
from pathlib import Path
from typing import Union
import zarr


def read_zarr(store: Union[str, Path], use_dask: bool = True) -> anndata.AnnData:
    """
    Assembles a reduced anndat instance based on a zarr store written with anndata.

    The anndata instance inclused .obs, .var, .uns and .X.
    In contrast to anndata._io.zarr.read_zarr(), this functon interfaces components of anndata as lazy dask arrays on
    the zarr store rather than loading the arrays from the zarr store fully into memory.

    :param store: Path to zarr group.
    :param use_dask: Whether to use lazy dask arrays where appropriate.
    :return: Anndata with components as dask array.
    """
    from anndata._io.zarr import _clean_uns, read_dataframe, read_attribute
    if isinstance(store, Path):
        store = str(store)

    f = zarr.open(store, mode="r")
    d = {
        "obs": read_dataframe(f["obs"]),  # TODO: could use a dask dataframe here.
        "var": read_dataframe(f["var"]),
        "uns": read_attribute(f["uns"]),
    }
    _clean_uns(d)
    adata = anndata.AnnData(**d)
    if use_dask:
        x = dask.array.from_zarr(url=store, component="X")
    else:
        x = read_attribute(f["X"])
    adata.X = x
    return adata
