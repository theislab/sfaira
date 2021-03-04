from typing import Union

from sfaira.data import DatasetGroup
from sfaira.data.dataloaders.super_group import DatasetSuperGroupSfaira


class DatasetGroupPancreas(DatasetGroup):

    def __init__(
        self,
        data_path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(data_path=data_path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "mouse_pancreas_2019_10x_pisco_001_10.1101/661728",
            "mouse_pancreas_2019_smartseq2_pisco_001_10.1101/661728",
            "mouse_pancreas_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001",
            "mouse_pancreas_2019_10x_thompson_001",
            "mouse_pancreas_2019_10x_thompson_002",
            "mouse_pancreas_2019_10x_thompson_003",
            "mouse_pancreas_2019_10x_thompson_004",
            "mouse_pancreas_2019_10x_thompson_005",
            "mouse_pancreas_2019_10x_thompson_006",
            "mouse_pancreas_2019_10x_thompson_007",
            "mouse_pancreas_2019_10x_thompson_008",
        ])
        super().__init__(datasets=dsg.flatten().datasets)
