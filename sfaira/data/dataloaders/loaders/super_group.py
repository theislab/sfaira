import pydoc
import os
from typing import Union

from sfaira.data import DatasetSuperGroup


class DatasetSuperGroupLoaders(DatasetSuperGroup):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
    ):
        """
        Class that sits ontop of a directory of data set directories that each contain a data set group.

        :param file_base:
        :param dir_prefix: Prefix to subselect directories by. Set to "" for no constraints.
        :param path:
        :param meta_path:
        :param cache_path:
        """
        # Directory choice hyperparamters:
        dir_prefix = "d"
        dir_exlcude = []
        # Collect all data loaders from files in directory:
        dataset_groups = []
        cwd = os.path.dirname(__file__)
        for f in os.listdir(cwd):
            if os.path.isdir(os.path.join(cwd, f)):  # only directories
                # Narrow down to data set directories:
                if f[:len(dir_prefix)] == dir_prefix and f not in dir_exlcude:
                    DatasetGroupDirectoryOriented = pydoc.locate(
                        "sfaira.sfaira.data.dataloaders.loaders." + f + ".DatasetGroupDirectoryOriented")
                    if DatasetGroupDirectoryOriented is not None:
                        dataset_groups.append(
                            DatasetGroupDirectoryOriented(path=path, meta_path=meta_path, cache_path=cache_path)
                        )
                    else:
                        print(f"WARNING: DatasetGroupDirectoryOriented was None for {f}")
        super().__init__(dataset_groups=dataset_groups)
