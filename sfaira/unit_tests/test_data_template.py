import os
import unittest

from sfaira.data import DatasetGroupDirectoryOriented


class TestDatasetTemplate(unittest.TestCase):
    dir_template: str = "./template_data"

    def test_load(self):
        """
        Address ToDos before running test to customize to your data set.
        :return:
        """
        remove_gene_version = True
        match_to_reference = None
        # ToDo: add correct module here as "YOUR_STUDY":
        # Addition coming soon: This path can either be in sfaira or in sfaira_extensions.
        # So far, this still has to be in sfaira.
        from sfaira.data.dataloaders.loaders.d10_1016_j_cmet_2019_01_021 import FILE_PATH
        ds = DatasetGroupDirectoryOriented(
            file_base=FILE_PATH,
            path=self.dir_template,
            meta_path=self.dir_template,
            cache_path=self.dir_template
        )
        # Test raw loading and caching:
        # You can set load_raw to True while debugging when caching works already to speed the test up,
        # but be sure to set load_raw to True for final tests.
        ds.load(
            remove_gene_version=remove_gene_version,
            match_to_reference=match_to_reference,
            load_raw=False,  # tests raw loading
            allow_caching=True  # tests caching
        )
        # Create cell type conversion table:
        for k, v in ds.datasets.items():
            v.load()
            # Write this directly into sfaira installation so that it can be committed via git.
            v.write_ontology_class_map(
                fn=os.path.join("/".join(FILE_PATH.split("/")[:-1]), v.fn_ontology_class_map_csv),
                protected_writing=False
            )
            # ToDo: conflicts are not automatically resolved, please go back to https://www.ebi.ac.uk/ols/ontologies/cl
            # for every mismatch or conflict and add the correct cell ontology class name into the .csv "target" column.
        # Test loading from cache:
        ds.load(
            remove_gene_version=remove_gene_version,
            match_to_reference=match_to_reference,
            load_raw=False,
            allow_caching=False
        )
        # Test concatenation:
        _ = ds.adata


if __name__ == '__main__':
    unittest.main()
