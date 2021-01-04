from typing import Union

from .external import DatasetGroupBase

from sfaira.data.dataloaders.mouse.d10_1101_661728 import Dataset as Dataset0001
from sfaira.data.dataloaders.mouse.d10_1101_661728 import Dataset as Dataset0002
from sfaira.data.dataloaders.mouse import Dataset as Dataset0003
from sfaira.data.dataloaders.mouse.d10_1016_j_cmet_2019_01_021.mouse_pancreas_2019_10x_thompson_001 import Dataset as Dataset0004
from sfaira.data.dataloaders.mouse import Dataset as Dataset0005
from sfaira.data.dataloaders.mouse import Dataset as Dataset0006
from sfaira.data.dataloaders.mouse import Dataset as Dataset0007
from sfaira.data.dataloaders.mouse import Dataset as Dataset0008
from sfaira.data.dataloaders.mouse.d10_1016_j_cmet_2019_01_021.mouse_pancreas_2019_10x_thompson_006 import Dataset as Dataset0009
from sfaira.data.dataloaders.mouse import Dataset as Dataset0010
from sfaira.data.dataloaders.mouse import Dataset as Dataset0011


class DatasetGroupPancreas(DatasetGroupBase):

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
            Dataset0007(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0008(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0009(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0010(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0011(path=path, meta_path=meta_path, cache_path=cache_path)
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.mouse import DatasetGroupPancreas
            self.datasets.update(DatasetGroupPancreas(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
