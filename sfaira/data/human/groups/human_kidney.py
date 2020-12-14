from typing import Union

from .external import DatasetGroupBase

from sfaira.data.human.d10_1038_s41467_019_10861_2.human_kidney_2019_10xSn_lake_001 import Dataset as Dataset0001
from sfaira.data.human.d10_1126_science_aat5031.human_kidney_2019_10x_stewart_001 import Dataset as Dataset0002
from sfaira.data.human.d10_1038_s41597_019_0351_8.human_kidney_2020_10x_liao_001 import Dataset as Dataset0003
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_kidney_2020_microwell_han_001 import Dataset as Dataset0004
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_kidney_2020_microwell_han_002 import Dataset as Dataset0005
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_kidney_2020_microwell_han_003 import Dataset as Dataset0006
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_kidney_2020_microwell_han_004 import Dataset as Dataset0007
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_kidney_2020_microwell_han_005 import Dataset as Dataset0008
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_kidney_2020_microwell_han_006 import Dataset as Dataset0009
from sfaira.data.human.d10_1038_s41586_020_2157_4.human_kidney_2020_microwell_han_007 import Dataset as Dataset0010


class DatasetGroupKidney(DatasetGroupBase):

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
            from sfaira_extension.data.human import DatasetGroupKidney
            self.datasets.update(DatasetGroupKidney(path=path, meta_path=meta_path).datasets)
        except ImportError:
            pass
