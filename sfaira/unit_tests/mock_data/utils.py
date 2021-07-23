import anndata
import scipy.sparse
import numpy as np
import os
import pandas as pd
from typing import Union

from sfaira.consts import AdataIdsSfaira, OCS
from sfaira.versions.genomes import GenomeContainer
from sfaira.versions.metadata import OntologyOboCustom

from .consts import ASSEMBLY_HUMAN, ASSEMBLY_MOUSE
from .loaders import DatasetSuperGroupMock

dir_store = os.path.join(os.path.dirname(__file__), "store")


def simulate_anndata(genes, n_obs, targets=None, assays=None, obo: Union[None, OntologyOboCustom] = None) -> \
        anndata.AnnData:
    """
    Simulate basic data example.

    :return: AnnData instance.
    """
    adata_ids_sfaira = AdataIdsSfaira()
    data = anndata.AnnData(
        np.random.randint(low=0, high=100, size=(n_obs, len(genes))).astype(np.float32)
    )
    if assays is not None:
        data.obs[adata_ids_sfaira.assay_sc] = [
            assays[np.random.randint(0, len(assays))]
            for _ in range(n_obs)
        ]
    if targets is not None:
        data.obs[adata_ids_sfaira.cellontology_class] = [
            targets[np.random.randint(0, len(targets))]
            for _ in range(n_obs)
        ]
        if obo is None:
            data.obs[adata_ids_sfaira.cellontology_id] = [
                OCS.cellontology_class.convert_to_id(x)
                if x not in [adata_ids_sfaira.unknown_celltype_identifier,
                             adata_ids_sfaira.not_a_cell_celltype_identifier]
                else x
                for x in data.obs[adata_ids_sfaira.cellontology_class].values
            ]
        else:
            data.obs[adata_ids_sfaira.cellontology_id] = [
                obo.convert_to_id(x)
                if x not in [adata_ids_sfaira.unknown_celltype_identifier,
                             adata_ids_sfaira.not_a_cell_celltype_identifier]
                else x
                for x in data.obs[adata_ids_sfaira.cellontology_class].values
            ]
    data.var[adata_ids_sfaira.gene_id_ensembl] = genes
    return data


def create_adata(celltypes, ncells, ngenes, assembly) -> anndata.AnnData:
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


def prepare_dsg() -> DatasetSuperGroupMock:
    """
    Prepares data set super group of mock data and return instance.
    """
    dsg = DatasetSuperGroupMock()
    return dsg


def prepare_store(store_format: str, rewrite: bool = False) -> str:
    """
    Prepares mock data store and returns path to store.
    """
    dir_store_formatted = dir_store + "_" + store_format
    adata_ids_sfaira = AdataIdsSfaira()
    dsg = prepare_dsg()
    if not os.path.exists(dir_store_formatted):
        os.mkdir(dir_store_formatted)
    # Only load files that are not already in cache.
    if rewrite:
        anticipated_files = np.unique([
            v.doi_journal for _, v in dsg.datasets.items()
        ]).tolist()
    else:
        anticipated_files = np.unique([
            v.doi_journal for _, v in dsg.datasets.items()
            if (not os.path.exists(os.path.join(dir_store_formatted, v.doi_cleaned_id + ".h5ad"))
                and store_format == "h5ad") or
               (not os.path.exists(os.path.join(dir_store_formatted, v.doi_cleaned_id)) and store_format == "dao")
        ]).tolist()
    dsg.subset(key=adata_ids_sfaira.doi_journal, values=anticipated_files)
    match_to_reference = {"human": ASSEMBLY_HUMAN, "mouse": ASSEMBLY_MOUSE}
    for ds in dsg.datasets.values():
        ds.load(allow_caching=True)
        ds.streamline_features(remove_gene_version=True, match_to_reference=match_to_reference,
                               subset_genes_to_type="protein_coding")
        ds.streamline_metadata(schema="sfaira", clean_obs=True, clean_var=True, clean_uns=True, clean_obs_names=True)
        if store_format == "dao":
            compression_kwargs = {"compressor": "default", "overwrite": True, "order": "C"}
        else:
            compression_kwargs = {}
        ds.write_distributed_store(dir_cache=dir_store_formatted, store_format=store_format, dense=True,
                                   chunks=128, compression_kwargs=compression_kwargs)
    return dir_store_formatted
