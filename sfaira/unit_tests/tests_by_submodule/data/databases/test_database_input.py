import os
import pytest

import dask.array
import numpy as np
import pandas as pd

from sfaira.consts import AdataIdsSfaira
from sfaira.data.dataloaders.databases.cellxgene import DatasetSuperGroupCellxgene
from sfaira.data.store.io.io_dao import read_dao
from sfaira.versions.genomes import GenomeContainer
from sfaira.unit_tests.data_for_tests.databases.utils import prepare_dsg_database
from sfaira.unit_tests.data_for_tests.databases.consts import CELLXGENE_DATASET_ID
from sfaira.unit_tests.data_for_tests.loaders import RELEASE_HUMAN, RELEASE_MOUSE
from sfaira.unit_tests.directories import DIR_DATA_DATABASES_CACHE, DIR_DATABASE_STORE_DAO

MATCH_TO_RELEASE = {"Homo sapiens": RELEASE_HUMAN, "Mus musculus": RELEASE_MOUSE}


@pytest.mark.parametrize("database", [
    ("cellxgene", ["id", CELLXGENE_DATASET_ID]),
])
@pytest.mark.parametrize("subset_genes_to_type", [None, "protein_coding"])
def test_streamline_features(database: str, subset_genes_to_type: str):
    """Check if feature streamlining from cellxgene is successful."""
    adata_ids_sfaira = AdataIdsSfaira()
    database, subset_args = database
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    assert len(dsg.dataset_groups) > 0
    dsg.load(clean_ontology_labels=True)
    organism = [v.organism for v in dsg.datasets.values()][0]
    gc = GenomeContainer(organism=organism, release=MATCH_TO_RELEASE[organism])
    # Define set of raw IDs:
    original_ids = dict([
        (k, np.array([x for x in v.adata.var.index.tolist() if x in gc.ensembl]))
        for k, v in dsg.datasets.items()])
    dsg.streamline_var(match_to_release=MATCH_TO_RELEASE,
                       remove_gene_version=True,
                       subset_genes_to_type=subset_genes_to_type)
    # Initialise reference gc to check target space inside of ds.
    if subset_genes_to_type is not None:
        gc.set(biotype=subset_genes_to_type)
    for k, v in dsg.datasets.items():
        if subset_genes_to_type is None:
            # Should have maintained original IDs.
            assert len(v.adata.var.index) == len(original_ids[k])
            assert np.all(v.adata.var.index == original_ids[k])
        else:
            # Should have expanded features to target space.
            assert np.all(v.adata.var[adata_ids_sfaira.feature_id].values == gc.ensembl)


@pytest.mark.parametrize("database", [
    ("cellxgene", ["id", CELLXGENE_DATASET_ID]),
])
@pytest.mark.parametrize("format", ["sfaira", ])
def test_streamline_metadata(database: str, format: str):
    """Check if meta data streamlining from cellxgene is successful."""
    database, subset_args = database
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load(clean_ontology_labels=True)
    dsg.streamline_var(match_to_release=MATCH_TO_RELEASE,
                       schema=format,
                       subset_genes_to_type="protein_coding")
    dsg.streamline_obs_uns(schema=format)
    adata = dsg.datasets[subset_args[1]].adata
    ids = AdataIdsSfaira()
    assert "CL:0000128" in adata.obs[ids.cell_type + ids.onto_id_suffix].values
    assert "oligodendrocyte" in adata.obs[ids.cell_type].values
    assert "MmusDv:0000061" in adata.obs[ids.development_stage + ids.onto_id_suffix].values
    assert "early adult stage" in adata.obs[ids.development_stage].values
    assert "UBERON:0002436" in adata.obs[ids.organ + ids.onto_id_suffix].values
    assert "primary visual cortex" in adata.obs[ids.organ].values


@pytest.mark.parametrize("store", ["dao", ])
@pytest.mark.parametrize("database", [
    ("cellxgene", ["id", CELLXGENE_DATASET_ID]),
])
def test_output_to_store(store: str, database: str):
    """Tests if store from cellxgene has correct content."""
    adata_ids_sfaira = AdataIdsSfaira()
    database, subset_args = database
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load(clean_ontology_labels=True)
    adata_original = dsg.datasets[CELLXGENE_DATASET_ID].adata.copy()
    dsg.streamline_var(match_to_release=MATCH_TO_RELEASE, schema="sfaira", subset_genes_to_type="protein_coding",
                       clean_var=True, subset_layer="counts")
    dsg.streamline_obs_uns(schema="sfaira", clean_obs=True, clean_uns=True, clean_obs_names=True,
                           keep_id_obs=True, keep_orginal_obs=False, keep_symbol_obs=True)
    dsg.write_distributed_store(dir_cache=DIR_DATABASE_STORE_DAO, store_format=store, dense=True)
    fn_store = os.path.join(DIR_DATABASE_STORE_DAO, subset_args[1])
    adata = read_dao(store=fn_store)
    ids = AdataIdsSfaira()
    assert "CL:0000128" in adata.obs[ids.cell_type + ids.onto_id_suffix].values
    assert "oligodendrocyte" in adata.obs[ids.cell_type].values
    assert "MmusDv:0000061" in adata.obs[ids.development_stage + ids.onto_id_suffix].values
    assert "early adult stage" in adata.obs[ids.development_stage].values
    # Check obs setup:
    anticipated_obs = adata_ids_sfaira.obs_keys
    anticipated_obs = [getattr(adata_ids_sfaira, x) for x in anticipated_obs]
    for k in adata_ids_sfaira.ontology_constrained:
        anticipated_obs.append(getattr(adata_ids_sfaira, k) + adata_ids_sfaira.onto_id_suffix)
    assert len(adata.obs.columns) == len(anticipated_obs)
    assert np.all([x in adata.obs.columns for x in anticipated_obs]), (adata.obs.columns, anticipated_obs)
    # Check var setup:
    organism = [v.organism for v in dsg.datasets.values()][0]
    gc = GenomeContainer(organism=organism, release=MATCH_TO_RELEASE[organism])
    gc.set(biotype="protein_coding")
    anticipated_var = adata_ids_sfaira.var_keys
    anticipated_var = [getattr(adata_ids_sfaira, x) for x in anticipated_var]
    assert len(adata.var.columns) == len(anticipated_var)
    assert np.all([x in adata.var.columns for x in anticipated_var]), (adata.var.columns, anticipated_var)
    assert adata.n_vars == gc.n_var
    assert np.all(adata.var_names == gc.ensembl)
    # Check X:
    # Compare a non-zero gene from .raw with entry written into store:
    # Note that adata::raw is written into store::X because subset_layer="counts" is used above.
    var_check = adata_original.raw.var_names[np.where(np.asarray(adata_original.raw.X.sum(axis=0)).flatten() > 0)[0]]
    var_check = [x for x in var_check if x in gc.ensembl][0]
    xj_original = np.asarray(adata_original.raw[:, var_check].X.todense()).flatten()
    x = dask.array.from_zarr(url=os.path.join(fn_store, "zarr"), component="X")
    var = pd.read_parquet(os.path.join(fn_store, "parquet", "var.parquet"), engine="pyarrow")
    xj = x[:, np.where(var["ensembl"].values == var_check)[0]].compute().flatten()
    assert np.sum(xj) > 0
    assert np.all(xj_original == xj)


def test_cellxgene_single_cell_subset():
    """
    This test can be used to debug cellxgene collection sub-setting.
    """
    dsg = DatasetSuperGroupCellxgene(data_path=DIR_DATA_DATABASES_CACHE, cache_metadata=True)
    assays = np.unique([str(x.assay_sc) for x in dsg.flatten().datasets.values() if x.assay_sc is not None])
    print(f"assay None: {np.sum([x.ncells for x in dsg.flatten().datasets.values() if x.assay_sc is None])} cells")
    for k in assays:
        print(f"assay {k}: {np.sum([x.ncells for x in dsg.flatten().datasets.values() if str(x.assay_sc) == k])} cells")
    n_total_cells = dsg.ncells
    dsg.subset(key="assay_sc", values="sfaira single cell library construction")
    print(f"found {n_total_cells} total cells and {dsg.ncells} selected cells")
    dsg.subset(key="assay_sc", values=["RNA assay", "10x technology", "sci-RNA-seq", "microwell-seq"])
    print(f"found {n_total_cells} total cells and {dsg.ncells} selected cells")
    dsg.subset(key="primary_data", values=True)
    print(f"found {n_total_cells} total cells and {dsg.ncells} selected cells")
    assert dsg.ncells > 1e7
