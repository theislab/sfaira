from typing import Union

from .external import DatasetGroup

from sfaira.data.dataloaders.loaders.d10_15252_embj_2018100811.human_eye_2019_10x_lukowski_001 import Dataset as Dataset0001
from sfaira.data.dataloaders.loaders.d10_1038_s41467_019_12780_8.human_eye_2019_10x_menon_001 import Dataset as Dataset0002
from sfaira.data.dataloaders.loaders.d10_1073_pnas_1914143116.human_eye_2019_10x_voigt_001 import Dataset as Dataset0003
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_eye_2020_microwell_han_001 import Dataset as Dataset0004


class DatasetGroupEye(DatasetGroup):

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
            from sfaira_extension.data.human import DatasetGroupEye
            self.datasets.update(DatasetGroupEye(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
