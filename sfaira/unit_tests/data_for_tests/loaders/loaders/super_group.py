import pydoc
import os
from typing import List
from warnings import warn
from sfaira.data import DatasetBase, DatasetGroup, DatasetSuperGroup, DatasetGroupDirectoryOriented

from sfaira.unit_tests.directories import DIR_DATA_LOADERS_CACHE


class DatasetSuperGroupMock(DatasetSuperGroup):
    """
    This is a DatasetSuperGroup which wraps the mock data loaders in the same directory.

    This class is designed to facilitate testing of code that requires data loaders without requiring raw data
    downloads as all mock data loaders operate on data that is simulated in the `load()` functions.
    A cache directory is established under ../cache.

    This class is a reduced and merged version of the sfaira loader super group class and the sfaira loader adapated
    DatasetGroupDirectoryOriented.
    """

    dataset_groups: List[DatasetGroupDirectoryOriented]

    def __init__(self):
        # Directory choice hyper-paramters:
        dir_prefix = "d"
        dir_exclude = []
        # Collect all data loaders from files in directory:
        dataset_groups = []
        cwd = os.path.dirname(__file__)
        for d in os.listdir(cwd):
            if os.path.isdir(os.path.join(cwd, d)):  # Iterate over mock studies (directories).
                if d[:len(dir_prefix)] == dir_prefix and d not in dir_exclude:  # Narrow down to data set directories
                    path_base = f"sfaira.unit_tests.data_for_tests.loaders.loaders.{d}"
                    path_dsg = pydoc.locate(f"{path_base}.FILE_PATH")
                    path_module = os.path.join(cwd, d)
                    for f in os.listdir(os.path.join(cwd, d)):  # Iterate over loaders in mock study (file).
                        datasets = []
                        if f.split(".")[-1] == "py" and f not in ["__init__.py"]:
                            file_module = ".".join(f.split(".")[:-1])
                            if path_dsg is not None:
                                load_func = pydoc.locate(f"{path_base}.{file_module}.load")
                                fn_yaml = os.path.join(path_module, file_module + ".yaml")
                                x = DatasetBase(
                                    data_path=DIR_DATA_LOADERS_CACHE,
                                    meta_path=DIR_DATA_LOADERS_CACHE,
                                    cache_path=DIR_DATA_LOADERS_CACHE,
                                    load_func=load_func,
                                    dict_load_func_annotation=None,
                                    sample_fn=None,
                                    sample_fns=None,
                                    yaml_path=fn_yaml,
                                )
                                x.read_ontology_class_map(fn=os.path.join(path_module, file_module + ".tsv"))
                                datasets.append(x)
                            else:
                                warn(f"DatasetGroupDirectoryOriented was None for {f}")
                        dataset_groups.append(DatasetGroup(datasets=dict([(x.id, x) for x in datasets]),
                                                           collection_id=d))
        super().__init__(dataset_groups=dataset_groups)
