import pytest
from typing import List

from sfaira.consts import AdataIdsSfaira
from sfaira.unit_tests.data_for_tests.databases.utils import prepare_dsg_database
from sfaira.unit_tests.data_for_tests.databases.consts import CELLXGENE_DATASET_ID
from sfaira.unit_tests.data_for_tests.loaders import ASSEMBLY_HUMAN, ASSEMBLY_MOUSE


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
    ids = AdataIdsSfaira()
    assert "CL:0000128" in adata.obs[ids.cell_type + ids.onto_id_suffix].values
    assert "oligodendrocyte" in adata.obs[ids.cell_type].values
    assert "HsapDv:0000087" in adata.obs[ids.development_stage + ids.onto_id_suffix].values
    assert "human adult stage" in adata.obs[ids.development_stage].values
    assert "UBERON:0000956" in adata.obs[ids.organ + ids.onto_id_suffix].values
    assert "cerebral cortex" in adata.obs[ids.organ].values