import os
from typing import Union

from .external import DatasetGroupBase

from .human_thyroid_2020_microwell_han_001 import Dataset as Dataset0001
from .human_thyroid_2020_microwell_han_002 import Dataset as Dataset0002


class DatasetGroupThyroid(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None
    ):
        datasets = [
            Dataset0001(path=path, meta_path=meta_path),
            Dataset0002(path=path, meta_path=meta_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.human import DatasetGroupThyroid
            self.datasets.update(DatasetGroupThyroid().datasets)
        except ImportError:
            pass
