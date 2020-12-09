import os
from typing import Union

from .external import DatasetGroupBase

from .human_mixed_2019_10x_szabo_001 import Dataset as Dataset0001


class DatasetGroupMixed(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None
    ):
        datasets = [
            Dataset0001(path=path, meta_path=meta_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            import sfaira_extension
            self.datasets.update(sfaira_extension.data.human.DatasetGroupMixed().datasets)
        except ImportError:
            pass
