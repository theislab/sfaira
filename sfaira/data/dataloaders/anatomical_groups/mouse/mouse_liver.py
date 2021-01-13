import os
from typing import Union

from .external import DatasetGroup

from sfaira.data.dataloaders.loaders.d10_1101_661728.mouse_liver_2019_10x_pisco_001 import Dataset as Dataset0001
from sfaira.data.dataloaders.loaders.d10_1101_661728.mouse_liver_2019_smartseq2_pisco_001 import Dataset as Dataset0002
from sfaira.data.dataloaders.loaders.d10_1016_j_cell_2018_02_001.mouse_liver_2018_microwell_han_001 import Dataset as Dataset0003
from sfaira.data.dataloaders.loaders.d10_1016_j_cell_2018_02_001.mouse_liver_2018_microwell_han_002 import Dataset as Dataset0004


class DatasetGroupLiver(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        datasets = [
            Dataset0001(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0002(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0003(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0004(path=path, meta_path=meta_path, cache_path=cache_path)
        ]
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.mouse import DatasetGroupLiver
            self.datasets.update(DatasetGroupLiver(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
