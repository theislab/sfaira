from typing import Union

from .external import DatasetGroup

from sfaira.data.dataloaders.loaders.d10_1038_s41467_019_12464_3.human_mixed_2019_10x_szabo_001 import Dataset as Dataset0001


class DatasetGroupMixed(DatasetGroup):

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
            from sfaira_extension.data.human import DatasetGroupMixed
            self.datasets.update(DatasetGroupMixed(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
