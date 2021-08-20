import os
import pydoc
import shutil
import re
from typing import Union

from sfaira.data import DatasetGroupDirectoryOriented, DatasetGroup, DatasetBase
from sfaira.data.utils import read_yaml
from sfaira.consts.utils import clean_doi
from sfaira.commands.questionary import sfaira_questionary

try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None


class DataloaderAnnotater:

    def __init__(self):
        self.WD = os.path.dirname(__file__)
        self.file_path = None
        self.file_path_sfairae = None
        self.meta_path = None
        self.cache_path = None
        self.dir_loader = None
        self.dir_loader_sfairae = None
        self.package_source = None

    def annotate(self, path_loader: str, path_data: str, doi: Union[str, None]):
        """
        Annotates a provided dataloader.

        Moderate the suggestions made here: Choose the best fit cell ontology label for your cells.
        Sfaira uses multiple mechanisms of finding matches, depending on how the free text was generated, these might be
        differentially successful. The proposed IDs groups are separate by ":|||:" strings to give you a visual anchor
        when going through these lists. You need to delete all of these division strings and all labels in the second
        columns other than the best fit label. Do not change the first column,
        (Note that columns are separated by ",")
        You can also manually check maps here: https://www.ebi.ac.uk/ols/ontologies/cl
        """
        if not doi:
            doi = sfaira_questionary(function='text',
                                     question='DOI:',
                                     default='10.1000/j.journal.2021.01.001')
            while not re.match(r'\b10\.\d+/[\w.]+\b', doi):
                print('[bold red]The entered DOI is malformed!')  # noqa: W605
                doi = sfaira_questionary(function='text',
                                         question='DOI:',
                                         default='10.1000/j.journal.2021.01.001')
        doi_sfaira_repr = clean_doi(doi)
        self._setup_loader(doi_sfaira_repr)
        self._annotate(path_data, path_loader, doi, doi_sfaira_repr)

    def _setup_loader(self, doi_sfaira_repr: str):
        """
        Define the file names, loader paths and base paths of loader collections for sfaira and sfaira_extension
        """
        dir_loader_sfaira = "sfaira.data.dataloaders.loaders."
        file_path_sfaira = "/" + "/".join(pydoc.locate(dir_loader_sfaira + "FILE_PATH").split("/")[:-1])
        if sfairae is not None:
            dir_loader_sfairae = "sfaira_extension.data.dataloaders.loaders."
            file_path_sfairae = "/" + "/".join(pydoc.locate(dir_loader_sfairae + "FILE_PATH").split("/")[:-1])
        else:
            file_path_sfairae = None
        # Check if loader name is a directory either in sfaira or sfaira_extension loader collections:
        if doi_sfaira_repr in os.listdir(file_path_sfaira):
            dir_loader = dir_loader_sfaira + "." + doi_sfaira_repr
            package_source = "sfaira"
        elif doi_sfaira_repr in os.listdir(file_path_sfairae):
            dir_loader = dir_loader_sfairae + "." + doi_sfaira_repr
            package_source = "sfairae"
        else:
            raise ValueError("data loader not found in sfaira and also not in sfaira_extension")
        file_path = pydoc.locate(dir_loader + ".FILE_PATH")
        meta_path = None
        cache_path = None
        # Clear dataset cache
        shutil.rmtree(cache_path, ignore_errors=True)

        self.file_path = file_path
        self.file_path_sfairae = file_path_sfairae
        self.meta_path = meta_path
        self.cache_path = cache_path
        self.dir_loader = dir_loader
        self.dir_loader_sfairae = None if sfairae is None else dir_loader_sfairae
        self.package_source = package_source

    def _get_ds(self, test_data: str):
        ds = DatasetGroupDirectoryOriented(
            file_base=self.file_path,
            data_path=test_data,
            meta_path=None,
            cache_path=None
        )

        return ds

    def buffered_load(self, test_data: str, doi_sfaira_repr: str):
        if not os.path.exists(test_data):
            raise ValueError(f"test-data directory {test_data} does not exist.")
        if doi_sfaira_repr not in os.listdir(test_data):
            raise ValueError(f"did not find data folder named {doi_sfaira_repr} in test-data directory "
                             f"{test_data}, only found {os.listdir(test_data)}")
        ds = self._get_ds(test_data=test_data)
        ds.load(
            remove_gene_version=False,
            match_to_reference=None,
            load_raw=True,  # Force raw load so non confound future tests by data loader bugs in previous versions.
            allow_caching=False,
            verbose=3
        )
        assert len(ds.ids) > 0, f"no data sets loaded, make sure raw data is in {test_data}, "\
                                f"found {os.listdir(os.path.join(test_data, doi_sfaira_repr))}"
        return ds

    def _annotate(self, test_data: str, path: str, doi: str, doi_sfaira_repr: str):
        ds = self.buffered_load(test_data=test_data, doi_sfaira_repr=doi_sfaira_repr)
        # Create cell type conversion table:
        cwd = os.path.dirname(self.file_path)
        dataset_module = str(cwd.split("/")[-1])
        # Group data sets by file module:
        # Note that if we were not grouping the cell type map .tsv files by file module, we could directly call
        # write_ontology_class_map on the ds.
        tsvs_written = []
        for f in os.listdir(cwd):
            if os.path.isfile(os.path.join(cwd, f)):  # only files
                # Narrow down to data set files:
                if f.split(".")[-1] == "py" and f.split(".")[0] not in ["__init__", "base", "group"]:
                    file_module = ".".join(f.split(".")[:-1])

                    # I) Instantiate Data set group to get all IDs of data sets associated with this .py file.
                    # Note that all data sets in this directory are already loaded in ds, so we just need the IDs.
                    DatasetFound = pydoc.locate(self.dir_loader + "." + file_module + ".Dataset")
                    # Load objects from name space:
                    # - load(): Loading function that return anndata instance.
                    # - SAMPLE_FNS: File name list for DatasetBaseGroupLoadingManyFiles
                    load_func = pydoc.locate(self.dir_loader + "." + file_module + ".load")
                    load_func_annotation = pydoc.locate(self.dir_loader + "." + file_module + ".LOAD_ANNOTATION")
                    # Also check sfaira_extension for additional load_func_annotation:
                    if self.package_source != "sfairae" and sfairae is not None:
                        load_func_annotation_sfairae = pydoc.locate(self.dir_loader_sfairae + "." + dataset_module +
                                                                    "." + file_module + ".LOAD_ANNOTATION")
                        # LOAD_ANNOTATION is a dictionary so we can use update to extend it.
                        if load_func_annotation_sfairae is not None and load_func_annotation is not None:
                            load_func_annotation.update(load_func_annotation_sfairae)
                        elif load_func_annotation_sfairae is not None and load_func_annotation is None:
                            load_func_annotation = load_func_annotation_sfairae
                    sample_fns = pydoc.locate(self.dir_loader + "." + file_module + ".SAMPLE_FNS")
                    fn_yaml = os.path.join(cwd, file_module + ".yaml")
                    fn_yaml = fn_yaml if os.path.exists(fn_yaml) else None
                    # Check for sample_fns in yaml:
                    if fn_yaml is not None:
                        assert os.path.exists(fn_yaml), f"did not find yaml {fn_yaml}"
                        yaml_vals = read_yaml(fn=fn_yaml)
                        if sample_fns is None and yaml_vals["meta"]["sample_fns"] is not None:
                            sample_fns = yaml_vals["meta"]["sample_fns"]
                    if sample_fns is None:
                        sample_fns = [None]
                    # Here we distinguish between class that are already defined and those that are not.
                    # The latter case arises if meta data are defined in YAMLs and _load is given as a function.
                    if DatasetFound is None:
                        datasets_f = [
                            DatasetBase(
                                data_path=test_data,
                                meta_path=self.meta_path,
                                cache_path=self.cache_path,
                                load_func=load_func,
                                dict_load_func_annotation=load_func_annotation,
                                sample_fn=x,
                                sample_fns=sample_fns if sample_fns != [None] else None,
                                yaml_path=fn_yaml,
                            ) for x in sample_fns
                        ]
                    else:
                        datasets_f = [
                            DatasetFound(
                                data_path=test_data,
                                meta_path=self.meta_path,
                                cache_path=self.cache_path,
                                load_func=load_func,
                                load_func_annotation=load_func_annotation,
                                sample_fn=x,
                                sample_fns=sample_fns if sample_fns != [None] else None,
                                yaml_path=fn_yaml,
                            ) for x in sample_fns
                        ]
                    # II) Build a data set group from the already loaded data sets and use the group ontology writing
                    # function.
                    dsg_f = DatasetGroup(datasets=dict([(x.id, ds.datasets[x.id]) for x in datasets_f]))
                    # III) Write this directly into the sfaira clone so that it can be committed via git.
                    # TODO any errors not to be caught here?
                    doi_sfaira_repr = f'd{doi.translate({ord(c): "_" for c in r"!@#$%^&*()[]/{};:,.<>?|`~-=_+"})}'
                    fn_tsv = os.path.join(path, doi_sfaira_repr, f"{file_module}.tsv")
                    dsg_f.write_ontology_class_map(
                        fn=fn_tsv,
                        protected_writing=True,
                        n_suggest=4,
                    )
                    tsvs_written.append(fn_tsv)
        print(f"Completed annotation. Wrote {len(tsvs_written)} files:\n" + "\n".join(tsvs_written))
