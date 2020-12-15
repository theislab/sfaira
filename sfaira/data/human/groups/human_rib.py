from typing import Union

from .external import DatasetGroupBase

from sfaira.data.human.d10_1038_s41586_020_2157_4.human_rib_2020_microwell_han_001 import Dataset as Dataset0001
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_rib_2020_microwell_han_002 import Dataset as Dataset0002


class DatasetGroupRib(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        super().__init__()
        datasets = [
            Dataset0001(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0002(path=path, meta_path=meta_path, cache_path=cache_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.human import DatasetGroupRib
            self.datasets.update(DatasetGroupRib(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
