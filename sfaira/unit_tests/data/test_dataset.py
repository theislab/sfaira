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


def test_dsgs_subset():
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.id
            assert v.organ == "lung", v.id


def test_dsgs_config_write_load():
    fn = dir_data + "/config.csv"
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.write_config(fn=fn)
    ds2 = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds2.load_config(fn=fn)
    assert np.all(ds.ids == ds2.ids)


def test_dsg_load():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load_all()


def test_dsg_adata():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    _ = ds.adata


"""
TODO tests from here on down require cached data for mouse lung
"""


def test_dsgs_adata():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    _ = ds.adata


def test_dsgs_load():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load_all()


@pytest.mark.parametrize("format", ["sfaira", "cellxgene"])
@pytest.mark.parametrize("clean", [True, False])
def test_dsgs_streamline(format: str, clean: bool):
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load_all()
    ds.streamline(format=format, clean=clean)


def test_dsg_load_backed_dense(genome="Mus_musculus_GRCm38_97"):
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load_all_tobacked(
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
    ds.load_all_tobacked(
        fn_backed=os.path.join(dir_data, 'test_backed_data.h5ad'),
        genome=genome,
        shuffled=False,
        as_dense=False,
        annotated_only=False
    )
    assert isinstance(ds.adata.X[:], scipy.sparse.csr_matrix), "%s" % type(ds.adata.X)
