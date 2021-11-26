import os

import anndata
import pytest

from sfaira.unit_tests.data_for_tests.loaders import RELEASE_HUMAN, PrepareData
from sfaira.unit_tests.directories import DIR_TEMP


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
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    if out_format == "cellxgene":
        # Other data data sets do not have complete enough annotation
        ds.subset(key="doi_journal", values=["no_doi_mock1", "no_doi_mock3"])
    ds.load()
    ds.streamline_features(remove_gene_version=False, match_to_release=RELEASE_HUMAN,
                           subset_genes_to_type=None)
    ds.streamline_metadata(schema=out_format, clean_obs=clean_obs, clean_var=clean_var,
                           clean_uns=clean_uns, clean_obs_names=clean_obs_names,
                           keep_id_obs=keep_id_obs, keep_orginal_obs=keep_orginal_obs, keep_symbol_obs=keep_symbol_obs)


@pytest.mark.parametrize("schema_version", ["2_0_0"])
@pytest.mark.parametrize("organism", ["Homo sapiens", "Mus musculus"])
def test_cellxgene_export(schema_version: str, organism: str):
    """

    This test can be extended by future versions.
    """
    from cellxgene_schema import validate

    class ValidatorInMemory(validate.Validator):
        """
        Helper class to validate adata in memory and raise errors as in error stream rather than outstream.

        The switch in log stream allows this test to be used as a unit test.
        """

        def validate_adata_inmemory(self, adata: anndata.AnnData):
            self.errors = []
            self.adata = adata
            self._set_schema_def()
            if not self.errors:
                self._deep_check()
            if self.warnings:
                self.warnings = ["WARNING: " + i for i in self.warnings]
            if self.errors:
                self.errors = ["ERROR: " + i for i in self.errors]
            if self.warnings or self.errors:
                print(self.warnings[:20])
                print(self.errors[:20])
                assert False

    ds = PrepareData().prepare_dsg(load=False)
    if organism == "Homo sapiens":
        ds.subset(key="doi_journal", values=["no_doi_mock1"])
    else:
        ds.subset(key="doi_journal", values=["no_doi_mock2"])
    ds.load()
    ds.streamline_features(schema="cellxgene:" + schema_version)
    ds.streamline_metadata(schema="cellxgene:" + schema_version, clean_obs=False, clean_var=True,
                           clean_uns=True, clean_obs_names=False,
                           keep_id_obs=True, keep_orginal_obs=False, keep_symbol_obs=True)
    counter = 0
    for ds in ds.datasets.values():
        # Validate in memory:
        # This directly tests the raw object.
        val = ValidatorInMemory()
        val.validate_adata_inmemory(adata=ds.adata)
        # Validate on disk:
        # This guards against potential issues that arise from h5 formatting.
        fn_temp = os.path.join(DIR_TEMP, "temp.h5ad")
        val = ValidatorInMemory()
        ds.adata.write_h5ad(filename=fn_temp)
        val.validate_adata(h5ad_path=fn_temp)
        counter += 1
    assert counter > 0, "no data sets to test"
