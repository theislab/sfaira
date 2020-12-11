import os
from typing import Union

from .external import DatasetGroupBase

from .mouse_kidney_2019_10x_pisco_001 import Dataset as Dataset0001
from .mouse_kidney_2019_smartseq2_pisco_001 import Dataset as Dataset0002
from .mouse_kidney_2018_microwell_han_001 import Dataset as Dataset0003
from .mouse_kidney_2018_microwell_han_002 import Dataset as Dataset0004


class DatasetGroupKidney(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None
    ):
        super().__init__()
        datasets = [
            Dataset0001(path=path, meta_path=meta_path),
            Dataset0002(path=path, meta_path=meta_path),
            Dataset0003(path=path, meta_path=meta_path),
            Dataset0004(path=path, meta_path=meta_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.mouse import DatasetGroupKidney
            self.datasets.update(DatasetGroupKidney(path=path, meta_path=meta_path).datasets)
        except ImportError:
            pass
