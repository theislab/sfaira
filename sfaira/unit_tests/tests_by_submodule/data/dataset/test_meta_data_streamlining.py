from cellxgene_schema.validate import validate_adata
import pytest

from sfaira.unit_tests.data_for_tests.loaders import ASSEMBLY_HUMAN, ASSEMBLY_MOUSE, prepare_dsg


@pytest.mark.parametrize("out_format", ["sfaira", "cellxgene"])
@pytest.mark.parametrize("clean_obs", [True, False])
@pytest.mark.parametrize("clean_var", [True, False])
@pytest.mark.parametrize("clean_uns", [True, False])
@pytest.mark.parametrize("clean_obs_names", [True, False])
@pytest.mark.parametrize("keep_id_obs", [True])
@pytest.mark.parametrize("keep_orginal_obs", [False])
@pytest.mark.parametrize("keep_symbol_obs", [True])
def test_dsgs_streamline_metadata(out_format: str, clean_obs: bool, clean_var: bool, clean_uns: bool,
                                  clean_obs_names: bool, keep_id_obs: bool, keep_orginal_obs: bool,
                                  keep_symbol_obs: bool):
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["human"])
    ds.subset(key="organ", values=["lung"])
    if out_format == "cellxgene":
        # Other data data sets do not have complete enough annotation
        ds.subset(key="doi_journal", values=["no_doi_mock1", "no_doi_mock3"])
    ds.load()
    ds.streamline_features(remove_gene_version=False, match_to_reference=ASSEMBLY_MOUSE,
                           subset_genes_to_type=None)
    ds.streamline_metadata(schema=out_format, clean_obs=clean_obs, clean_var=clean_var,
                           clean_uns=clean_uns, clean_obs_names=clean_obs_names,
                           keep_id_obs=keep_id_obs, keep_orginal_obs=keep_orginal_obs, keep_symbol_obs=keep_symbol_obs)


@pytest.mark.parametrize("schema_version", ["1_1_0"])
@pytest.mark.parametrize("organism", ["human", "mouse"])
def test_cellxgene_export(schema_version: str, organism: str):
    """

    This test can be extended by future versions.
    """
    ds = prepare_dsg(load=False)
    if organism == "human":
        ds.subset(key="doi_journal", values=["no_doi_mock1"])
    else:
        ds.subset(key="doi_journal", values=["no_doi_mock2"])
    ds.load()
    ds.streamline_features(remove_gene_version=False,
                           match_to_reference={"human": ASSEMBLY_HUMAN, "mouse": ASSEMBLY_MOUSE},
                           subset_genes_to_type=None)
    ds.streamline_metadata(schema="cellxgene:" + schema_version, clean_obs=False, clean_var=True,
                           clean_uns=True, clean_obs_names=False,
                           keep_id_obs=True, keep_orginal_obs=False, keep_symbol_obs=True)
    counter = 0
    for ds in ds.datasets.values():
        validate_adata(adata=ds.adata, shallow=False)
        counter += 1
    assert counter > 0, "no data sets to test"
