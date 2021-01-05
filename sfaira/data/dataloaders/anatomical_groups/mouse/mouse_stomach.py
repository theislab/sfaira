import os
from typing import Union

from .external import DatasetGroup

from sfaira.data.dataloaders.loaders.d10_1016_j_cell_2018_02_001.mouse_stomach_2018_microwell_han_001 import Dataset as Dataset0001


class DatasetGroupStomach(DatasetGroup):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        datasets = [
            Dataset0001(path=path, meta_path=meta_path, cache_path=cache_path)
        ]
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.mouse import DatasetGroupStomach
            self.datasets.update(DatasetGroupStomach(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
