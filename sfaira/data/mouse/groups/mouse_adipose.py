import os
from typing import Union

from .external import DatasetGroupBase

from sfaira.data.mouse.d10_1101_661728.mouse_adipose_2019_10x_pisco_001 import Dataset as Dataset0001
from sfaira.data.mouse.d10_1101_661728.mouse_adipose_2019_smartseq2_pisco_001 import Dataset as Dataset0002
from sfaira.data.mouse.d10_1101_661728.mouse_adipose_2019_smartseq2_pisco_002 import Dataset as Dataset0003
from sfaira.data.mouse.d10_1101_661728.mouse_adipose_2019_smartseq2_pisco_003 import Dataset as Dataset0004
from sfaira.data.mouse.d10_1101_661728.mouse_adipose_2019_smartseq2_pisco_004 import Dataset as Dataset0005


class DatasetGroupAdipose(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        super().__init__()
        datasets = [
            Dataset0001(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0002(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0003(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0004(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0005(path=path, meta_path=meta_path, cache_path=cache_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.mouse import DatasetGroupAdipose
            self.datasets.update(DatasetGroupAdipose(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
