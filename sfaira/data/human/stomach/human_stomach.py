import os
from typing import Union

from .external import DatasetGroupBase

from .human_stomach_2020_microwell_han_001 import Dataset as Dataset0001
from .human_stomach_2020_microwell_han_002 import Dataset as Dataset0002
from .human_stomach_2020_microwell_han_003 import Dataset as Dataset0003
from .human_stomach_2020_microwell_han_004 import Dataset as Dataset0004
from .human_stomach_2020_microwell_han_005 import Dataset as Dataset0005
from .human_stomach_2020_microwell_han_006 import Dataset as Dataset0006
from .human_stomach_2020_microwell_han_007 import Dataset as Dataset0007
from .human_stomach_2020_microwell_han_008 import Dataset as Dataset0008
from .human_stomach_2020_microwell_han_009 import Dataset as Dataset0009
from .human_stomach_2020_microwell_han_010 import Dataset as Dataset0010


class DatasetGroupStomach(DatasetGroupBase):

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
            Dataset0004(path=path, meta_path=meta_path),
            Dataset0005(path=path, meta_path=meta_path),
            Dataset0006(path=path, meta_path=meta_path),
            Dataset0007(path=path, meta_path=meta_path),
            Dataset0008(path=path, meta_path=meta_path),
            Dataset0009(path=path, meta_path=meta_path),
            Dataset0010(path=path, meta_path=meta_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.human import DatasetGroupStomach
            self.datasets.update(DatasetGroupStomach(path=path, meta_path=meta_path).datasets)
        except ImportError:
            pass
