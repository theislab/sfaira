import unittest

from sfaira.data import DatasetGroupDirectoryOriented


class TestDatasetTemplate(unittest.TestCase):
    dir_data: str = "./test_data"
    dir_meta: str = "./test_data/meta"

    def test_load(self):
        """
        Address ToDos before running test to customize to your data set.
        :return:
        """
        celltype_version = None
        remove_gene_version = True
        match_to_reference = None
        # ToDo: add correct module here as "YOUR_STUDY":
        from sfaira.data.dataloaders.loaders.YOUR_STUDY import FILE_PATH
        ds = DatasetGroupDirectoryOriented(
            file_base=FILE_PATH,
            path=self.dir_data,
            meta_path=self.dir_meta,
            cache_path=self.dir_data
        )
        # Test raw loading and caching:
        ds.load(
            celltype_version=celltype_version,
            fn=None,
            remove_gene_version=remove_gene_version,
            match_to_reference=match_to_reference,
            load_raw=True,  # tests raw loading
            allow_caching=True  # tests caching
        )
        # Test loading from cache:
        ds.load(
            celltype_version=celltype_version,
            fn=None,
            remove_gene_version=remove_gene_version,
            match_to_reference=match_to_reference,
            load_raw=False,
            allow_caching=False
        )
        # Test concatenation:
        _ = ds.adata


if __name__ == '__main__':
    unittest.main()
