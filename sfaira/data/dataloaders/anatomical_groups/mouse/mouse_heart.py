from typing import Union

from sfaira.data import DatasetGroup
from sfaira.data.dataloaders.super_group import Universe


class DatasetGroupHeart(DatasetGroup):

    def __init__(
        self,
        data_path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = Universe(data_path=data_path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "mouse_heart_2019_10x_pisco_001_10.1101/661728",
            "mouse_heart_2019_smartseq2_pisco_001_10.1101/661728",
            "mouse_heart_2019_smartseq2_pisco_002_10.1101/661728"
        ])
        super().__init__(datasets=dsg.flatten().datasets)
