import anndata
import scipy.sparse
import numpy as np
import os
import pandas as pd
import pathlib
from sfaira.versions.genomes import GenomeContainer

from sfaira.unit_tests.directories import DIR_DATA_LOADERS_CACHE, DIR_DATA_LOADERS_STORE_DAO, \
    DIR_DATA_LOADERS_STORE_H5AD, save_delete
from .consts import ASSEMBLY_HUMAN, ASSEMBLY_MOUSE
from .loaders import DatasetSuperGroupMock

MATCH_TO_REFERENCE = {"human": ASSEMBLY_HUMAN, "mouse": ASSEMBLY_MOUSE}


def _create_adata(celltypes, ncells, ngenes, assembly) -> anndata.AnnData:
    """
    Usesd by mock data loaders.
    """
    gc = GenomeContainer(assembly=assembly)
    gc.subset(biotype="protein_coding")
    genes = gc.ensembl[:ngenes]
    x = scipy.sparse.csc_matrix(np.random.randint(low=0, high=100, size=(ncells, ngenes)))
    var = pd.DataFrame(index=genes)
    obs = pd.DataFrame({}, index=["cell_" + str(i) for i in range(ncells)])
    if len(celltypes) > 0:
        obs["free_annotation"] = [celltypes[i] for i in np.random.choice(len(celltypes), size=ncells, replace=True)]
    adata = anndata.AnnData(X=x, obs=obs, var=var)
    return adata


def _load_script(dsg, rewrite: bool, match_to_reference):
    dsg.load(allow_caching=True, load_raw=rewrite)
    dsg.streamline_features(remove_gene_version=True, match_to_reference=match_to_reference)
    dsg.streamline_metadata(schema="sfaira", clean_obs=True, clean_var=True, clean_uns=True, clean_obs_names=True)
    return dsg


def prepare_dsg(rewrite: bool = False, load: bool = True) -> DatasetSuperGroupMock:
    """
    Prepares data set super group of mock data and returns instance.

    Use this do testing involving a data set group.
    """
    # Make sure cache exists:
    if not os.path.exists(DIR_DATA_LOADERS_CACHE):
        pathlib.Path(DIR_DATA_LOADERS_CACHE).mkdir(parents=True, exist_ok=True)
    dsg = DatasetSuperGroupMock()
    if load:
        dsg = _load_script(dsg=dsg, rewrite=rewrite, match_to_reference=MATCH_TO_REFERENCE)
    return dsg


def prepare_store(store_format: str, rewrite: bool = False, rewrite_store: bool = False) -> str:
    """
    Prepares mock data store and returns path to store.

    Use this do testing involving a data set store.
    """
    dir_store_formatted = {
        "dao": DIR_DATA_LOADERS_STORE_DAO,
        "h5ad": DIR_DATA_LOADERS_STORE_H5AD,
    }[store_format]
    if not os.path.exists(dir_store_formatted):
        pathlib.Path(dir_store_formatted).mkdir(parents=True, exist_ok=True)
    dsg = prepare_dsg(rewrite=rewrite, load=False)
    for k, ds in dsg.datasets.items():
        if store_format == "dao":
            compression_kwargs = {"compressor": "default", "overwrite": True, "order": "C"}
        else:
            compression_kwargs = {}
        if store_format == "dao":
            anticipated_fn = os.path.join(dir_store_formatted, ds.doi_cleaned_id)
        elif store_format == "h5ad":
            anticipated_fn = os.path.join(dir_store_formatted, ds.doi_cleaned_id + ".h5ad")
        else:
            assert False
        if rewrite_store and os.path.exists(anticipated_fn):
            # Can't write if h5ad already exists.
            # Delete store to writing if forced.
            save_delete(anticipated_fn)
        # Only rewrite if necessary
        if rewrite_store or not os.path.exists(anticipated_fn):
            ds = _load_script(dsg=ds, rewrite=rewrite, match_to_reference=MATCH_TO_REFERENCE)
            ds.write_distributed_store(dir_cache=dir_store_formatted, store_format=store_format, dense=True,
                                       chunks=128, compression_kwargs=compression_kwargs)
    return dir_store_formatted
