import os
from typing import Union

from .external import DatasetGroupBase

from .mouse_diaphragm_2019_smartseq2_pisco_001 import Dataset as Dataset0001


class DatasetGroupDiaphragm(DatasetGroupBase):

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
            from sfaira_extension.data.mouse import DatasetGroupDiaphragm
            self.datasets.update(DatasetGroupDiaphragm(path=path, meta_path=meta_path).datasets)
        except ImportError:
            pass
