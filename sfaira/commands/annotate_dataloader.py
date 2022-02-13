import numpy as np
import os
import pydoc
from rich import print
import shutil
from typing import Union

from sfaira.commands.utils import get_pydoc
from sfaira.consts.utils import clean_doi
from sfaira.data import DatasetGroupDirectoryOriented, DatasetGroup, DatasetBase
from sfaira.data.utils import read_yaml

try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None


class DataloaderAnnotater:

    def __init__(self):
        self.WD = os.path.dirname(__file__)
        self.file_path = None
        self.meta_path = None
        self.cache_path = None
        self.pydoc_handle = None

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
        doi_sfaira_repr = clean_doi(doi)
        self._setup_loader(path_loader, doi_sfaira_repr)
        self._annotate(path_data, path_loader, doi, doi_sfaira_repr)

    def _setup_loader(self, path_loader: str, doi_sfaira_repr: str):
        """
        Define the file names, loader paths and base paths of loader collections for sfaira and sfaira_extension
        """
        file_path, pydoc_handle = get_pydoc(path_loader=path_loader, doi_sfaira_repr=doi_sfaira_repr)
        meta_path = None
        cache_path = None
        # Clear dataset cache
        shutil.rmtree(cache_path, ignore_errors=True)

        self.file_path = file_path
        self.meta_path = meta_path
        self.cache_path = cache_path
        self.pydoc_handle = pydoc_handle

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
                    DatasetFound = pydoc.locate(self.pydoc_handle + "." + file_module + ".Dataset")
                    # Load objects from name space:
                    # - load(): Loading function that return anndata instance.
                    # - SAMPLE_FNS: File name list for DatasetBaseGroupLoadingManyFiles
                    load_func = pydoc.locate(self.pydoc_handle + "." + file_module + ".load")
                    load_func_annotation = pydoc.locate(self.pydoc_handle + "." + file_module + ".LOAD_ANNOTATION")
                    sample_fns = pydoc.locate(self.pydoc_handle + "." + file_module + ".SAMPLE_FNS")
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
                    doi_sfaira_repr = clean_doi(doi)
                    fn_tsv = os.path.join(path, doi_sfaira_repr, f"{file_module}")
                    # Define .tsvs to write:
                    attrs = [
                        k for k in dsg_f._adata_ids.ontology_constrained
                        if np.any([getattr(v, k + "_obs_key") is not None for v in dsg_f.datasets.values()])
                    ]
                    dsg_f.write_ontology_class_maps(
                        fn=fn_tsv,
                        attrs=attrs,
                        protected_writing=True,
                        n_suggest=4,
                    )
                    tsvs_written.append((fn_tsv, attrs))
        print("[bold blue]Completed annotation.")
        print('[bold orange]Sfaira butler: "Up next, follow these steps until the next call of the sfaira CLI:"')
        self.action_counter = 1
        print(f'[bold orange]               "{self.action_counter}) Proceed to chose ontology symbols for each free '
              f'text label in the tsv files:"')
        for prefix, attrs in tsvs_written:
            for attr in attrs:
                print(f'[bold orange]                    -{prefix}_{attr}.tsv')
        print('[bold orange]                "Each tsv has two columns: free text labels found in the data on the left '
              'and suggestions on the right."')
        print('[bold orange]                "Each suggested symbol lies between two : characters."')
        print('[bold orange]                ": is a separator between suggested symbols and :|||: between symbol '
              'groups that were found through different search strategies."')
        print('[bold orange]                "Take care to not remove the tab separators in the table."')
        print('[bold orange]                "You only need to finish the second column now."')
        print('[bold orange]                "The third column with ontology IDs is added by `finalize` in phase 3."')
        self.action_counter += 1
        print(f'[bold orange]               "{self.action_counter}) Then proceed to finish .yaml file if not already '
              f'done."')
        self.action_counter += 1
        print(f'[bold orange]               "{self.action_counter}) Then proceed to phase 3 '
              f'sfaira finalize-dataloader."')
        self.action_counter += 1
