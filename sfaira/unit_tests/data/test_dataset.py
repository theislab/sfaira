import numpy as np
import os
import scipy.sparse

from sfaira.data import DatasetSuperGroup
from sfaira.data import DatasetSuperGroupSfaira

dir_data = "../test_data"
dir_meta = "../test_data/meta"


def test_instantiate():
    _ = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)


def test_load():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load_all()


def test_adata():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["bladder"])
    _ = ds.adata


def test_load():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load_all()


def test_adata():
    ds = DatasetSuperGroupSfaira(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    _ = ds.adata


def test_load_backed_dense(genome="Mus_musculus_GRCm38_97"):
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


def test_load_backed_sparse(genome="Mus_musculus_GRCm38_97"):
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
