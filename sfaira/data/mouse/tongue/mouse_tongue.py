import os
from typing import Union

from .external import DatasetGroupBase

from .mouse_tongue_2019_10x_pisco_001 import Dataset as Dataset0001
from .mouse_tongue_2019_smartseq2_pisco_001 import Dataset as Dataset0002


class DatasetGroupTongue(DatasetGroupBase):

    def __init__(
        self, 
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None
    ):
        datasets = [
            Dataset0001(path=path, meta_path=meta_path),
            Dataset0002(path=path, meta_path=meta_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            import sfaira_extension.api as sfairae
            datasets.update(sfairae.data.mouse.DatasetGroupTongue().datasets)
        except ImportError:
            pass
