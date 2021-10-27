import pydoc
import os
from typing import List, Union
from warnings import warn
from sfaira.data.dataloaders.base.dataset_group import DatasetSuperGroup, DatasetGroupDirectoryOriented


class DatasetSuperGroupLoaders(DatasetSuperGroup):

    dataset_groups: List[DatasetGroupDirectoryOriented]

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
    ):
        """
        Class that sits on top of a directory of data set directories that each contain a data set group.

        :param file_base:
        :param dir_prefix: Prefix to sub-select directories by. Set to "" for no constraints.
        :param data_path:
        :param meta_path:
        :param cache_path:
        """
        # Directory choice hyperparamters:
        dir_prefix = "d"
        dir_exclude = []
        # Collect all data loaders from files in directory:
        dataset_groups = []
        cwd = os.path.dirname(__file__)
        for f in os.listdir(cwd):
            if os.path.isdir(os.path.join(cwd, f)):  # only directories
                if f[:len(dir_prefix)] == dir_prefix and f not in dir_exclude:  # Narrow down to data set directories
                    path_dsg = str(pydoc.locate(f"sfaira.data.dataloaders.loaders.{f}.FILE_PATH"))
                    if path_dsg is not None:
                        try:
                            dsg = DatasetGroupDirectoryOriented(
                                file_base=path_dsg,
                                data_path=data_path,
                                meta_path=meta_path,
                                cache_path=cache_path
                            )
                            dsg.collection_id = f
                            dataset_groups.append(dsg)
                        except IndexError as e:
                            raise IndexError(f"{e} for '{cwd}', '{f}', '{path_dsg}'")
                    else:
                        warn(f"DatasetGroupDirectoryOriented was None for {f}")
        super().__init__(dataset_groups=dataset_groups)
