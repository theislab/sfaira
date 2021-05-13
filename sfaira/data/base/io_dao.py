import anndata
import dask.array
import dask.dataframe
import os
import pandas as pd
from pathlib import Path
import pickle
from typing import List, Tuple, Union
import zarr


def _buffered_path(path_base, path, fn):
    path_base = os.path.join(path_base, path)
    if not os.path.exists(path_base):
        os.makedirs(path_base)
    return os.path.join(path_base, fn)


def path_obs(path):
    return _buffered_path(path_base=path, path="parquet", fn="obs.parquet")


def path_var(path):
    return _buffered_path(path_base=path, path="parquet", fn="var.parquet")


def path_uns(path):
    return _buffered_path(path_base=path, path="pickle", fn="uns.pickle")


def path_x(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return os.path.join(path, "zarr")


def write_dao(store: Union[str, Path], adata, chunks, compression_kwargs):
    """
    Writes a distributed access optimised ("dao") store of a dataset based on an AnnData instance.

    The following components are saved:
        - .X: as zarr array which can be interfaced with zarr or dask (or xarray).
        - .obs: as parquet table which can be interfaced with pandas.DataFrame (and dask.dataframe.DataFrame).
        - .var: as parquet table which can be interfaced with pandas.DataFrame (and dask.dataframe.DataFrame).
        - .uns: as a pickle to be flexible with values here.

    TODO: If layers become relevant for this store, they can be added into the zarr group.
    TODO: If obsp, varp become relevant for this store, they can be added into the zarr group.
    TODO: If obsm, varm become relevant for this store, they can be added into the zarr group.

    :param store: File name of the store (zarr group).
    :param adata: Anndata to save.
    :param chunks: Chunking of .X for zarr.
    :param compression_kwargs: Compression kwargs for zarr.
    """
    # Write numeric matrix as zarr array:
    f = zarr.open(path_x(store), mode="w")
    f.create_dataset("X", data=adata.X, chunks=chunks, **compression_kwargs)
    # Write .uns into pickle:
    with open(path_uns(store), "wb") as f:
        pickle.dump(obj=adata.uns, file=f)
    # Write .obs and .var as a separate file as this can be easily interfaced with DataFrames.
    adata.obs.to_parquet(path=path_obs(store), engine='pyarrow', compression='snappy', index=None)
    adata.var.to_parquet(path=path_var(store), engine='pyarrow', compression='snappy', index=None)


def read_dao(store: Union[str, Path], use_dask: bool = True, columns: Union[None, List[str]] = None,
             obs_separate: bool = False) -> \
        Union[Tuple[anndata.AnnData, Union[dask.dataframe.DataFrame, pd.DataFrame]], anndata.AnnData]:
    """
    Assembles an AnnData instance based on distributed access optimised ("dao") store of a dataset.

    See write_distributed_access_optimised() for the expected format of the store.
    In particular, if use_dask is True:

        - .X is interfaced as a dask Array

    Can return representation of .obs separately, which makes sense if a HPC framework is used for this tabular format
    which is not supported by AnnData as an .obs entry.

    :param store: Path to zarr group.
    :param use_dask: Whether to use lazy dask arrays where appropriate.
    :param columns: Which columns to read into the obs copy in the output, see pandas.read_parquet().
    :param obs_separate: Whether to return .obs as a separate return value or in the returned AnnData.
    :return: Tuple of:
        - AnnData with .X as dask array.
        - obs table separately as dataframe
    """
    if use_dask:
        x = dask.array.from_zarr(url=path_x(store), component="X")
    else:
        f = zarr.open(path_x(store), mode="r")
        x = f["X"]  # Select member of group.
        x = x[...]  # Load into memory.
    # Read pickle:
    with open(path_uns(store), "rb") as f:
        uns = pickle.load(file=f)
    # Read tables:
    obs = pd.read_parquet(path_obs(store), columns=columns, engine="pyarrow")
    var = pd.read_parquet(path_var(store), engine="pyarrow")
    # Convert to categorical variables where possible to save memory:
    for k, dtype in zip(list(obs.columns), obs.dtypes):
        if dtype == "object":
            obs[k] = obs[k].astype(dtype="category")
    d = {"var": var, "uns": uns}
    # Assemble AnnData without obs to save memory:
    adata = anndata.AnnData(**d, shape=x.shape)
    # Need to add these attributes after initialisation so that they are not evaluated:
    adata.X = x
    if obs_separate:
        return adata, obs
    else:
        adata.obs = obs
        return adata
