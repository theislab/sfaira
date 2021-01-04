from typing import Union

from .external import DatasetGroupBase

from sfaira.data.dataloaders.loaders.d10_1038_nmeth_4407.human_brain_2017_DroNcSeq_habib_001 import Dataset as Dataset0001
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_brain_2020_microwell_han_001 import Dataset as Dataset0002
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_brain_2020_microwell_han_002 import Dataset as Dataset0003
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_brain_2020_microwell_han_003 import Dataset as Dataset0004
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_brain_2020_microwell_han_004 import Dataset as Dataset0005
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_brain_2020_microwell_han_005 import Dataset as Dataset0006
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_brain_2020_microwell_han_006 import Dataset as Dataset0007


class DatasetGroupBrain(DatasetGroupBase):

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
            Dataset0005(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0006(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0007(path=path, meta_path=meta_path, cache_path=cache_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.human import DatasetGroupBrain
            self.datasets.update(DatasetGroupBrain(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
