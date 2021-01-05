from typing import Union

from .external import DatasetGroup

from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_epityphlon_2020_microwell_han_001 import Dataset as Dataset0001


class DatasetGroupEpityphlon(DatasetGroup):

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
            from sfaira_extension.data.human import DatasetGroupEpityphlon
            self.datasets.update(DatasetGroupEpityphlon(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
