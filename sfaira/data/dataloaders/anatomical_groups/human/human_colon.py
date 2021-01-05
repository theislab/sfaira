from typing import Union

from .external import DatasetGroup

from sfaira.data.dataloaders.loaders.d10_1016_j_cell_2018_08_067.human_colon_2019_10x_kinchen_001 import Dataset as Dataset0001
from sfaira.data.dataloaders.loaders.d10_1016_j_cell_2019_06_029.human_colon_2019_10x_smilie_001 import Dataset as Dataset0002
from sfaira.data.dataloaders.loaders.d10_1084_jem_20191130.human_colon_2019_10x_wang_001 import Dataset as Dataset0003
from sfaira.data.dataloaders.loaders.d10_1038_s41590_020_0602_z.human_colon_2020_10x_james_001 import Dataset as Dataset0004
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_colon_2020_microwell_han_001 import Dataset as Dataset0005
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_colon_2020_microwell_han_002 import Dataset as Dataset0006
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_colon_2020_microwell_han_003 import Dataset as Dataset0007
from sfaira.data.dataloaders.loaders.d10_1038_s41586_020_2157_4.human_colon_2020_microwell_han_004 import Dataset as Dataset0008


class DatasetGroupColon(DatasetGroup):

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
            Dataset0004(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0005(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0006(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0007(path=path, meta_path=meta_path, cache_path=cache_path),
            Dataset0008(path=path, meta_path=meta_path, cache_path=cache_path)
        ]
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
        # Load versions from extension if available:
        try:
            from sfaira_extension.data.human import DatasetGroupColon
            self.datasets.update(DatasetGroupColon(path=path, meta_path=meta_path, cache_path=cache_path).datasets)
        except ImportError:
            pass
