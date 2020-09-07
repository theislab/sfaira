import os
from typing import Union

from .external import DatasetGroupBase

from .human_spleen_2019_10x_madissoon_001 import Dataset as Dataset0001
from .human_spleen_2020_microwell_han_001 import Dataset as Dataset0002
from .human_spleen_2020_microwell_han_002 import Dataset as Dataset0003


class DatasetGroupSpleen(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None
    ):
        datasets = [
            Dataset0001(path=path, meta_path=meta_path),
            Dataset0002(path=path, meta_path=meta_path),
            Dataset0003(path=path, meta_path=meta_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            import sfaira_extension as sfairae
            datasets.update(sfairae.data.human.DatasetGroupSpleen().datasets)
        except ImportError:
            pass
