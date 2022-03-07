import os
import pytest

from sfaira.consts import AdataIdsSfaira
from sfaira.data.store.io.io_dao import read_dao
from sfaira.unit_tests.data_for_tests.databases.utils import prepare_dsg_database
from sfaira.unit_tests.data_for_tests.databases.consts import CELLXGENE_DATASET_ID
from sfaira.unit_tests.data_for_tests.loaders import RELEASE_HUMAN, RELEASE_MOUSE
from sfaira.unit_tests.directories import DIR_DATABASE_STORE_DAO

MATCH_TO_RELEASE = {"Homo sapiens": RELEASE_HUMAN, "Mus musculus": RELEASE_MOUSE}


@pytest.mark.parametrize("database", [
    ("cellxgene", ["id", CELLXGENE_DATASET_ID]),
])
@pytest.mark.parametrize("subset_genes_to_type", [None, "protein_coding", ])
def test_streamline_features(database: str, subset_genes_to_type: str):
    database, subset_args = database
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load()
    dsg.streamline_features(match_to_release=MATCH_TO_RELEASE, subset_genes_to_type=subset_genes_to_type)


@pytest.mark.parametrize("database", [
    ("cellxgene", ["id", CELLXGENE_DATASET_ID]),
])
@pytest.mark.parametrize("format", ["sfaira", ])
def test_streamline_metadata(database: str, format: str):
    database, subset_args = database
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load()
    dsg.streamline_features(match_to_release=MATCH_TO_RELEASE, subset_genes_to_type="protein_coding")
    dsg.streamline_metadata(schema=format)
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
    database, subset_args = database
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load()
    dsg.streamline_features(match_to_release=MATCH_TO_RELEASE, subset_genes_to_type="protein_coding")
    dsg.streamline_metadata(schema="sfaira", clean_obs=True, clean_uns=True, clean_var=True, clean_obs_names=True,
                            keep_id_obs=True, keep_orginal_obs=False, keep_symbol_obs=True)
    dsg.write_distributed_store(dir_cache=DIR_DATABASE_STORE_DAO, store_format=store, dense=True)
    fn_store = os.path.join(DIR_DATABASE_STORE_DAO, subset_args[1])
    adata = read_dao(store=fn_store)
    ids = AdataIdsSfaira()
    assert "CL:0000128" in adata.obs[ids.cell_type + ids.onto_id_suffix].values
    assert "oligodendrocyte" in adata.obs[ids.cell_type].values
    assert "MmusDv:0000061" in adata.obs[ids.development_stage + ids.onto_id_suffix].values
    assert "early adult stage" in adata.obs[ids.development_stage].values
