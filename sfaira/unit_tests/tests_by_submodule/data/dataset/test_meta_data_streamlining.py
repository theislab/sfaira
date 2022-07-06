import os

import anndata
import numpy as np
import pytest

from sfaira.consts.adata_fields import AdataIdsSfaira, AdataIdsCellxgene_v2_0_0
from sfaira.unit_tests.data_for_tests.loaders import RELEASE_HUMAN, PrepareData, PrepareDataExport
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
    if out_format == "cellxgene":
        ds = PrepareDataExport().prepare_dsg(load=False)
    else:
        ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_var(schema=out_format, remove_gene_version=False, match_to_release=RELEASE_HUMAN,
                      subset_genes_to_type=None, clean_var=clean_var)
    ds.streamline_obs_uns(schema=out_format, clean_obs=clean_obs,
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
        k = "no_doi_mock1"
    else:
        k = "no_doi_mock2"
    ds.subset(key="doi_journal", values=[k])
    ds.load()
    adata_pre = ds.adata_ls[0].copy()
    ds.streamline_var(match_to_release=None, schema="cellxgene:" + schema_version, clean_var=True)
    ds.streamline_obs_uns(schema="cellxgene:" + schema_version, clean_obs=False,
                          clean_uns=True, clean_obs_names=False,
                          keep_id_obs=True, keep_orginal_obs=False, keep_symbol_obs=True)
    adata_post = ds.adata_ls[0].copy()
    # Custom checks:
    adata_ids_sfaira = AdataIdsSfaira()
    if schema_version == "2_0_0":
        adata_ids_cellxgene = AdataIdsCellxgene_v2_0_0()
    else:
        assert False
    # .obs
    # Check that number of observations and obs names were not changed:
    assert adata_pre.obs.shape[0] == adata_post.obs.shape[0]
    assert np.all(adata_post.obs_names == adata_pre.obs_names)
    # Check that original columns are removed:
    assert np.all([adata_ids_cellxgene.onto_original_suffix not in k for k in adata_post.obs.columns])
    assert np.all([adata_ids_sfaira.onto_original_suffix not in k for k in adata_post.obs.columns])
    # .var
    # Check that feature space remained unchanged:
    assert adata_pre.var.shape[0] == adata_post.var.shape[0]
    assert np.all(adata_post.var.index == adata_pre.var.index)
    # .X, layers
    # Check that sum over count matrix remained unchanged:
    assert np.asarray(np.sum(adata_pre.X)).flatten() == \
           np.asarray(np.sum(adata_post.raw.X)).flatten()
    # .uns
    # None yet.
    # Native cellxgene checks:
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
