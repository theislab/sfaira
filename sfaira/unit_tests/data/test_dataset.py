import numpy as np
import os
import pytest
import scipy.sparse

from sfaira.data import DatasetSuperGroup
from sfaira.data import DatasetSuperGroupSfaira

dir_data = "../test_data"
dir_meta = "../test_data/meta"


def test_dsgs_instantiate():
    _ = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)


@pytest.mark.parametrize("organ", ["intestine", "ileum"])
def test_dsgs_subset_dataset_wise(organ: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=[organ])
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.organism
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ


def test_dsgs_config_write_load():
    fn = dir_data + "/config.csv"
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.write_config(fn=fn)
    ds2 = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds2.load_config(fn=fn)
    assert np.all(ds.ids == ds2.ids)


"""
TODO tests from here on down require cached data for mouse lung
"""


def test_dsgs_adata():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load(remove_gene_version=True)
    _ = ds.adata


def test_dsgs_load():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load(remove_gene_version=False)


@pytest.mark.parametrize("organ", ["lung"])
@pytest.mark.parametrize("celltype", ["T cell"])
def test_dsgs_subset_cell_wise(organ: str, celltype: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=[organ])
    ds.load(remove_gene_version=False)
    ds.subset_cells(key="cellontology_class", values=celltype)
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.id
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ
            for y in np.unique(v.adata.obs[v._adata_ids_sfaira.cell_ontology_class].values):
                assert v.ontology_container_sfaira.cellontology_class.is_a(query=y, reference=celltype), y


@pytest.mark.parametrize("out_format", ["sfaira", "cellxgene"])
@pytest.mark.parametrize("clean_obs", [True, False])
@pytest.mark.parametrize("clean_var", [True, False])
@pytest.mark.parametrize("clean_uns", [True, False])
def test_dsgs_streamline(out_format: str, clean_obs: bool, clean_var: bool, clean_uns: bool):
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load(remove_gene_version=True)
    ds.streamline(format=out_format, allow_uns_sfaira=False, clean_obs=clean_obs, clean_var=clean_var, clean_uns=clean_uns)


def test_dsg_load_backed_dense(genome="Mus_musculus_GRCm38_97"):
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load_tobacked(
        fn_backed=os.path.join(dir_data, 'test_backed_data.h5ad'),
        genome=genome,
        shuffled=True,
        as_dense=True,
        annotated_only=False
    )
    assert isinstance(ds.adata.X[:], np.ndarray), "%s" % type(ds.adata.X)


def test_dsg_load_backed_sparse(genome="Mus_musculus_GRCm38_97"):
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load_tobacked(
        fn_backed=os.path.join(dir_data, 'test_backed_data.h5ad'),
        genome=genome,
        shuffled=False,
        as_dense=False,
        annotated_only=False
    )
    assert isinstance(ds.adata.X[:], scipy.sparse.csr_matrix), "%s" % type(ds.adata.X)


def test_dsg_load():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load()


def test_dsg_adata():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    _ = ds.adata
