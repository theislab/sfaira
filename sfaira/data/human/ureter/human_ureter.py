import os
from typing import Union

from .external import DatasetGroupBase

from .human_ureter_2020_microwell_han_001 import Dataset as Dataset0001


class DatasetGroupUreter(DatasetGroupBase):

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
            from sfaira_extension.data.human import DatasetGroupUreter
            self.datasets.update(DatasetGroupUreter(path=path, meta_path=meta_path).datasets)
        except ImportError:
            pass
