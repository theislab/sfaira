import numpy as np
import os
import scipy.sparse
import unittest

from sfaira.data import mouse, DatasetSuperGroup


class TestDatasetGroups(unittest.TestCase):
    dir_data: str = "./test_data"
    dir_meta: str = "./test_data/meta"

    def test_load(self):
        ds = mouse.DatasetGroupLung(path=self.dir_data, meta_path=self.dir_meta)
        ds.load_all()

    def test_adata(self):
        ds = mouse.DatasetGroupBladder(path=self.dir_data, meta_path=self.dir_meta)
        _ = ds.adata


class TestDatasetSuperGroups(unittest.TestCase):
    dir_data: str = "./test_data"
    dir_meta: str = "./test_data/meta"

    def test_load(self):
        ds = DatasetSuperGroup(
            dataset_groups=[
                mouse.DatasetGroupLung(path=self.dir_data, meta_path=self.dir_meta)
            ]
        )
        ds.load_all()

    def test_adata(self):
        ds = DatasetSuperGroup(
            dataset_groups=[
                mouse.DatasetGroupLung(path=self.dir_data, meta_path=self.dir_meta)
            ]
        )
        _ = ds.adata

    def test_load_backed_dense(self, genome="Mus_musculus_GRCm38_97"):
        ds = DatasetSuperGroup(
            dataset_groups=[
                mouse.DatasetGroupLung(path=self.dir_data, meta_path=self.dir_meta)
            ]
        )
        ds.load_all_tobacked(
            fn_backed=os.path.join(self.dir_data, 'test_backed_data.h5ad'),
            genome=genome,
            shuffled=True,
            as_dense=True,
            annotated_only=False
        )
        assert isinstance(ds.adata.X[:], np.ndarray), "%s" % type(ds.adata.X)

    def test_load_backed_sparse(self, genome="Mus_musculus_GRCm38_97"):
        ds = DatasetSuperGroup(
            dataset_groups=[
                mouse.DatasetGroupLung(path=self.dir_data, meta_path=self.dir_meta)
            ]
        )
        ds.load_all_tobacked(
            fn_backed=os.path.join(self.dir_data, 'test_backed_data.h5ad'),
            genome=genome,
            shuffled=False,
            as_dense=False,
            annotated_only=False
        )
        assert isinstance(ds.adata.X[:], scipy.sparse.csr_matrix), "%s" % type(ds.adata.X)


if __name__ == '__main__':
    unittest.main()
