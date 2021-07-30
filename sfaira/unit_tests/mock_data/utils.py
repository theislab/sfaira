import anndata
import scipy.sparse
import numpy as np
import os
import pandas as pd
from sfaira.versions.genomes import GenomeContainer
from .consts import ASSEMBLY_HUMAN, ASSEMBLY_MOUSE
from .loaders import DatasetSuperGroupMock

dir_store = os.path.join(os.path.dirname(__file__), "store")


def _create_adata(celltypes, ncells, ngenes, assembly) -> anndata.AnnData:
    """
    Usesd by mock data loaders.
    """
    gc = GenomeContainer(assembly=assembly)
    gc.subset(biotype="protein_coding")
    genes = gc.ensembl[:ngenes]
    x = scipy.sparse.csc_matrix(np.random.randint(low=0, high=100, size=(ncells, ngenes)))
    var = pd.DataFrame(index=genes)
    obs = pd.DataFrame({
        "free_annotation": [celltypes[i] for i in np.random.choice(a=[0, 1], size=ncells, replace=True)]
    }, index=["cell_" + str(i) for i in range(ncells)])
    adata = anndata.AnnData(X=x, obs=obs, var=var)
    return adata


def prepare_dsg(rewrite: bool = False) -> DatasetSuperGroupMock:
    """
    Prepares data set super group of mock data and return instance.

    Use this do testing involving a data set group.
    """
    dsg = DatasetSuperGroupMock()
    match_to_reference = {"human": ASSEMBLY_HUMAN, "mouse": ASSEMBLY_MOUSE}
    dsg.load(allow_caching=True, load_raw=rewrite)
    dsg.streamline_features(remove_gene_version=True, match_to_reference=match_to_reference)
    dsg.streamline_metadata(schema="sfaira", clean_obs=True, clean_var=True, clean_uns=True, clean_obs_names=True)
    return dsg


def prepare_store(store_format: str, rewrite: bool = False) -> str:
    """
    Prepares mock data store and returns path to store.

    Use this do testing involving a data set store.
    """
    dir_store_formatted = dir_store + "_" + store_format
    dsg = prepare_dsg(rewrite=rewrite)
    if not os.path.exists(dir_store_formatted):
        os.mkdir(dir_store_formatted)
    for ds in dsg.datasets.values():
        if store_format == "dao":
            compression_kwargs = {"compressor": "default", "overwrite": True, "order": "C"}
        else:
            compression_kwargs = {}
        ds.write_distributed_store(dir_cache=dir_store_formatted, store_format=store_format, dense=True,
                                   chunks=128, compression_kwargs=compression_kwargs)
    return dir_store_formatted
