import numpy as np
import os
import pytest

from sfaira.data import DatasetSuperGroup
from sfaira.data import Universe

from ..mock_data import ASSEMBLY_MOUSE, DIR_TEMP, prepare_dsg


def test_dsgs_instantiate():
    _ = Universe(data_path=".", meta_path=".", cache_path=".")


@pytest.mark.parametrize("organ", ["intestine", "ileum"])
def test_dsgs_subset_dataset_wise(organ: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=[organ])
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.organism
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ


def test_dsgs_config_write_load():
    fn = os.path.join(DIR_TEMP + "config.csv")
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.write_config(fn=fn)
    ds2 = prepare_dsg()
    ds2.load_config(fn=fn)
    assert np.all(ds.ids == ds2.ids)


"""
TODO tests from here on down require cached data for mouse lung
"""


def test_dsgs_adata():
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    _ = ds.adata


def test_dsgs_load():
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()


@pytest.mark.parametrize("organ", ["lung"])
@pytest.mark.parametrize("celltype", ["T cell"])
def test_dsgs_subset_cell_wise(organ: str, celltype: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=[organ])
    ds.load()
    ds.subset_cells(key="cellontology_class", values=celltype)
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.id
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ
            for y in np.unique(v.adata.obs[v._adata_ids.cellontology_class].values):
                assert v.ontology_container_sfaira.cellontology_class.is_a(query=y, reference=celltype), y


@pytest.mark.parametrize("out_format", ["sfaira", "cellxgene"])
@pytest.mark.parametrize("uns_to_obs", [True, False])
@pytest.mark.parametrize("clean_obs", [True, False])
@pytest.mark.parametrize("clean_var", [True, False])
@pytest.mark.parametrize("clean_uns", [True, False])
@pytest.mark.parametrize("clean_obs_names", [True, False])
def test_dsgs_streamline_metadata(out_format: str, uns_to_obs: bool, clean_obs: bool, clean_var: bool, clean_uns: bool,
                                  clean_obs_names: bool):
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=False, match_to_reference=ASSEMBLY_MOUSE,
                           subset_genes_to_type=None)
    ds.streamline_metadata(schema=out_format, clean_obs=clean_obs, clean_var=clean_var,
                           clean_uns=clean_uns, clean_obs_names=clean_obs_names)


@pytest.mark.parametrize("match_to_reference", ["Mus_musculus.GRCm38.102", {"mouse": ASSEMBLY_MOUSE}])
@pytest.mark.parametrize("remove_gene_version", [False, True])
@pytest.mark.parametrize("subset_genes_to_type", [None, "protein_coding"])
def test_dsgs_streamline_features(match_to_reference: str, remove_gene_version: bool, subset_genes_to_type: str):
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=remove_gene_version, match_to_reference=match_to_reference,
                           subset_genes_to_type=subset_genes_to_type)


@pytest.mark.parametrize("store", ["h5ad"])
@pytest.mark.parametrize("dense", [False])
@pytest.mark.parametrize("clean_obs", [False, True])
def test_dsg_write_store(store: str, dense: bool, clean_obs: bool):
    _ = prepare_dsg()


def test_dsg_load():
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load()


def test_dsg_adata():
    ds = prepare_dsg()
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    _ = ds.adata
