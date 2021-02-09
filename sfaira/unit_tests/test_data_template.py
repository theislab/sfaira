import os
import pydoc
import unittest

from sfaira.data import DatasetGroupDirectoryOriented, DatasetGroup


class TestDatasetTemplate(unittest.TestCase):

    dir_template: str = "./template_data"

    def test_load(self):
        """
        Unit test to assist with data set contribution.

        The workflow for contributing a data set with this data loader is as follows:

        1. Write a data loader and add it into the loader directory of your local sfaira installation.
        2. Address ToDos below.
        3. Run this unit test until you are not getting errors from your data loader anymore.

        In the process of this unit test, this data loader will have written putative cell type maps from your
        annotation to the cell ontology.

        4. Moderate the suggestions made here: Choose the best fit cell ontology label for your cells.
        Sfaira uses multiple mechanisms of finding matches, depending on how the free text was generated, these might be
        differentially successfull. The proposed IDs groups are separate by ":|||:" strings to give you a visial anchor
        when going through these lists. You need to delete all of these division strings and all labels in the second
        columns other than the best fit label. Do not change the first column,
        (Note that columns are separated by ",")
        You can also manually check maps here: https://www.ebi.ac.uk/ols/ontologies/cl
        5. Run this unit test for a last time to check the cell type maps.

        :return:
        """
        remove_gene_version = True
        match_to_reference = None
        classmap_by_file = True  # ToDo build one class map per file or per data loader (potentially many per file)
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
        cwd = os.path.dirname(FILE_PATH)
        dataset_module = str(cwd.split("/")[-1])
        if classmap_by_file:
            for f in os.listdir(cwd):
                if os.path.isfile(os.path.join(cwd, f)):  # only files
                    # Narrow down to data set files:
                    if f.split(".")[-1] == "py" and f.split(".")[0] not in ["__init__", "base", "group"]:
                        file_module = ".".join(f.split(".")[:-1])
                        DatasetFound = pydoc.locate(
                            "sfaira.data.dataloaders.loaders." + dataset_module + "." +
                            file_module + ".Dataset")
                        # Check if global objects are available:
                        # - SAMPLE_FNS: for DatasetBaseGroupLoadingManyFiles
                        # - SAMPLE_IDS: for DatasetBaseGroupLoadingOneFile
                        sample_fns = pydoc.locate(
                            "sfaira.data.dataloaders.loaders." + dataset_module + "." +
                            file_module + ".SAMPLE_FNS")
                        sample_ids = pydoc.locate(
                            "sfaira.data.dataloaders.loaders." + dataset_module + "." +
                            file_module + ".SAMPLE_IDS")
                        if sample_fns is not None and sample_ids is None:
                            # DatasetBaseGroupLoadingManyFiles:
                            datasets_f = [
                                DatasetFound(
                                    sample_fn=x,
                                    path=self.dir_template,
                                    meta_path=self.dir_template,
                                    cache_path=self.dir_template
                                )
                                for x in sample_fns
                            ]
                        elif sample_fns is None and sample_ids is not None:
                            # DatasetBaseGroupLoadingManyFiles:
                            datasets_f = [
                                DatasetFound(
                                    sample_id=x,
                                    path=self.dir_template,
                                    meta_path=self.dir_template,
                                    cache_path=self.dir_template
                                )
                                for x in sample_ids
                            ]
                        elif sample_fns is not None and sample_ids is not None:
                            raise ValueError(f"sample_fns and sample_ids both found for {f}")
                        else:
                            datasets_f = [DatasetFound(
                                path=self.dir_template,
                                meta_path=self.dir_template,
                                cache_path=self.dir_template
                            )]
                        # Build a data set group from the already loaded data sets and use the group ontology writing
                        # function.
                        current_ids = [x.id for x in datasets_f]
                        dsg_f = DatasetGroup(datasets=dict([(x, ds.datasets[x]) for x in current_ids]))
                        # Write this directly into sfaira installation so that it can be committed via git.
                        dsg_f.write_ontology_class_map(
                            fn=os.path.join(cwd, file_module + ".csv"),
                            protected_writing=True,
                            n_suggest=4,
                        )
        else:
            for k, v in ds.datasets.items():
                # Write this directly into sfaira installation so that it can be committed via git.
                v.write_ontology_class_map(
                    fn=os.path.join("/".join(FILE_PATH.split("/")[:-1]), v.fn_ontology_class_map_csv),
                    protected_writing=True,
                    n_suggest=10,
                )

        # ToDo: conflicts are not automatically resolved, please go back to
        #  https://www.ebi.ac.uk/ols/ontologies/cl
        #  for every mismatch or conflict and add the correct cell ontology class name into the .csv "target" column.

        # Test loading from cache:
        ds = DatasetGroupDirectoryOriented(
            file_base=FILE_PATH,
            path=self.dir_template,
            meta_path=self.dir_template,
            cache_path=self.dir_template
        )
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
