import os
from typing import Union

from .external import DatasetGroupBase

from .human_brain_2017_DroNcSeq_habib_001 import Dataset as Dataset0001
from .human_brain_2020_microwell_han_001 import Dataset as Dataset0002
from .human_brain_2020_microwell_han_002 import Dataset as Dataset0003
from .human_brain_2020_microwell_han_003 import Dataset as Dataset0004
from .human_brain_2020_microwell_han_004 import Dataset as Dataset0005
from .human_brain_2020_microwell_han_005 import Dataset as Dataset0006
from .human_brain_2020_microwell_han_006 import Dataset as Dataset0007


class DatasetGroupBrain(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None
    ):
        datasets = [
            Dataset0001(path=path, meta_path=meta_path),
            Dataset0002(path=path, meta_path=meta_path),
            Dataset0003(path=path, meta_path=meta_path),
            Dataset0004(path=path, meta_path=meta_path),
            Dataset0005(path=path, meta_path=meta_path),
            Dataset0006(path=path, meta_path=meta_path),
            Dataset0007(path=path, meta_path=meta_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            import sfaira_extension.api as sfairae
            datasets.update(sfairae.data.human.DatasetGroupBrain().datasets)
        except ImportError:
            pass
