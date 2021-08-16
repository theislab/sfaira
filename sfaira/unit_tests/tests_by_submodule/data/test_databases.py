import os
import pytest
import shutil
from typing import List

from sfaira.consts import AdataIdsCellxgene
from sfaira.unit_tests.directories import DIR_DATA_DATABASES_CACHE
from sfaira.unit_tests.data_for_tests.databases.utils import prepare_dsg_database
from sfaira.unit_tests.data_for_tests.databases.consts import CELLXGENE_DATASET_ID
from sfaira.unit_tests.data_for_tests.loaders import ASSEMBLY_HUMAN, ASSEMBLY_MOUSE


# Execute this one first so that data sets are only downloaded once. Named test_a for this reason.
@pytest.mark.parametrize("database", ["cellxgene", ])
@pytest.mark.parametrize("subset_args", [None, ["id", CELLXGENE_DATASET_ID], ])
def test_a_dsgs_download(database: str, subset_args: List[str]):
    """
    Tests if downloading of data base entries works.

    Warning, deletes entire database unit test cache.
    """
    if os.path.exists(DIR_DATA_DATABASES_CACHE):
        shutil.rmtree(DIR_DATA_DATABASES_CACHE)
    dsg = prepare_dsg_database(database=database, download=False)
    if subset_args is not None:
        dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.download()


@pytest.mark.parametrize("database", ["cellxgene", ])
@pytest.mark.parametrize("subset_args", [["id", CELLXGENE_DATASET_ID], ["organism", "human"], ])
def test_dsgs_subset(database: str, subset_args: List[str]):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])


@pytest.mark.parametrize("database", ["cellxgene", ])
@pytest.mark.parametrize("subset_args", [["id", CELLXGENE_DATASET_ID], ])
@pytest.mark.parametrize("match_to_reference", [None, {"human": ASSEMBLY_HUMAN, "mouse": ASSEMBLY_MOUSE}, ])
def test_dsgs_adata(database: str, subset_args: List[str], match_to_reference: dict):
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load()
    if match_to_reference is not None:
        dsg.streamline_features(remove_gene_version=True, match_to_reference=match_to_reference)
        dsg.streamline_metadata(schema="sfaira", clean_obs=True, clean_var=True, clean_uns=True, clean_obs_names=True)
    _ = dsg.adata


@pytest.mark.parametrize("database", ["cellxgene", ])
@pytest.mark.parametrize("subset_args", [["id", CELLXGENE_DATASET_ID], ])
@pytest.mark.parametrize("match_to_reference", [{"human": ASSEMBLY_HUMAN, "mouse": ASSEMBLY_MOUSE}, ])
@pytest.mark.parametrize("subset_genes_to_type", [None, "protein_coding", ])
def test_dsgs_streamline_features(database: str, subset_args: List[str], match_to_reference: dict,
                                  subset_genes_to_type: str):
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load()
    dsg.streamline_features(match_to_reference=match_to_reference, subset_genes_to_type=subset_genes_to_type)


@pytest.mark.parametrize("database", ["cellxgene", ])
@pytest.mark.parametrize("subset_args", [["id", CELLXGENE_DATASET_ID], ])
@pytest.mark.parametrize("format", ["sfaira", ])
def test_dsgs_streamline_metadata(database: str, subset_args: List[str], format: str):
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.load()
    dsg.streamline_features(match_to_reference={"human": ASSEMBLY_HUMAN, "mouse": ASSEMBLY_MOUSE},
                            subset_genes_to_type="protein_coding")
    dsg.streamline_metadata(schema=format)
    adata = dsg.datasets[subset_args[1]].adata
    ids = AdataIdsCellxgene()
    assert "CL:0000128" in adata.obs[ids.cell_type + ids.ontology_id_suffix].values
    assert "oligodendrocyte" in adata.obs[ids.cell_type].values
    assert "HsapDv:0000087" in adata.obs[ids.development_stage + ids.ontology_id_suffix].values
    assert "human adult stage" in adata.obs[ids.development_stage].values
    assert "UBERON:0000956" in adata.obs[ids.organ + ids.ontology_id_suffix].values
    assert "cerebral cortex" in adata.obs[ids.organ].values
