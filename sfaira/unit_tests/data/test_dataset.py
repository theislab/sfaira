import numpy as np
import os
import pytest

from sfaira.data import DatasetSuperGroup
from sfaira.data import Universe

MOUSE_GENOME_ANNOTATION = "Mus_musculus.GRCm38.102"

dir_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data")
dir_meta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data/meta")


def test_dsgs_instantiate():
    _ = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)


@pytest.mark.parametrize("organ", ["intestine", "ileum"])
def test_dsgs_subset_dataset_wise(organ: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=[organ])
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.organism
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ


def test_dsgs_config_write_load():
    fn = dir_data + "/config.csv"
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.write_config(fn=fn)
    ds2 = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds2.load_config(fn=fn)
    assert np.all(ds.ids == ds2.ids)


"""
TODO tests from here on down require cached data for mouse lung
"""


def test_dsgs_adata():
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    _ = ds.adata


def test_dsgs_load():
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()


@pytest.mark.parametrize("organ", ["lung"])
@pytest.mark.parametrize("celltype", ["T cell"])
def test_dsgs_subset_cell_wise(organ: str, celltype: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
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
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=False, match_to_reference=MOUSE_GENOME_ANNOTATION,
                           subset_genes_to_type=None)
    ds.streamline_metadata(schema=out_format, uns_to_obs=uns_to_obs, clean_obs=clean_obs, clean_var=clean_var,
                           clean_uns=clean_uns, clean_obs_names=clean_obs_names)


@pytest.mark.parametrize("match_to_reference", ["Mus_musculus.GRCm38.102", {"mouse": MOUSE_GENOME_ANNOTATION}])
@pytest.mark.parametrize("remove_gene_version", [False, True])
@pytest.mark.parametrize("subset_genes_to_type", [None, "protein_coding"])
def test_dsgs_streamline_features(match_to_reference: str, remove_gene_version: bool, subset_genes_to_type: str):
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=remove_gene_version, match_to_reference=match_to_reference,
                           subset_genes_to_type=subset_genes_to_type)


@pytest.mark.parametrize("store", ["h5ad"])
@pytest.mark.parametrize("dense", [False])
@pytest.mark.parametrize("clean_obs", [False, True])
def test_dsg_write_store(store: str, dense: bool, clean_obs: bool):
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=True, match_to_reference={"mouse": MOUSE_GENOME_ANNOTATION},
                           subset_genes_to_type="protein_coding")
    ds.streamline_metadata(schema="sfaira", uns_to_obs=False, clean_obs=clean_obs, clean_var=True, clean_uns=True,
                           clean_obs_names=True)
    ds.write_distributed_store(dir_cache=os.path.join(dir_data, "store"), store_format=store, dense=dense)


def test_dsg_load():
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load()


def test_dsg_adata():
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    _ = ds.adata
