from typing import Union

from .external import DatasetGroupBase

from .human_lung_2019_10x_braga_001 import Dataset as Dataset0001
from .human_lung_2019_10x_braga_002 import Dataset as Dataset0002
from .human_lung_2019_10x_madissoon_001 import Dataset as Dataset0003
from .human_lung_2019_dropseq_braga_003 import Dataset as Dataset0004
from .human_lung_2020_10x_habermann_001 import Dataset as Dataset0005
from .human_lung_2020_10x_lukassen_001 import Dataset as Dataset0006
from .human_lung_2020_10x_lukassen_002 import Dataset as Dataset0007
from .human_lung_2020_10x_miller_001 import Dataset as Dataset0008
from .human_lung_2020_10x_travaglini_001 import Dataset as Dataset0009
from .human_lung_2020_microwell_han_001 import Dataset as Dataset0010
from .human_lung_2020_microwell_han_002 import Dataset as Dataset0011
from .human_lung_2020_microwell_han_003 import Dataset as Dataset0012
from .human_lung_2020_microwell_han_004 import Dataset as Dataset0013
from .human_lung_2020_microwell_han_005 import Dataset as Dataset0014
from .human_lung_2020_smartseq2_travaglini_002 import Dataset as Dataset0015


class DatasetGroupLung(DatasetGroupBase):

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
            Dataset0007(path=path, meta_path=meta_path),
            Dataset0008(path=path, meta_path=meta_path),
            Dataset0009(path=path, meta_path=meta_path),
            Dataset0010(path=path, meta_path=meta_path),
            Dataset0011(path=path, meta_path=meta_path),
            Dataset0012(path=path, meta_path=meta_path),
            Dataset0013(path=path, meta_path=meta_path),
            Dataset0014(path=path, meta_path=meta_path),
            Dataset0015(path=path, meta_path=meta_path),
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))
        # Load versions from extension if available:
        try:
            import sfaira_extension
            self.datasets.update(sfaira_extension.data.human.DatasetGroupLung().datasets)
        except ImportError:
            pass
