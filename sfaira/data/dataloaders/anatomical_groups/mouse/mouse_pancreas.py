from typing import Union

from sfaira.data import DatasetGroup, DatasetSuperGroupSfaira


class DatasetGroupPancreas(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(path=path, meta_path=meta_path, cache_path=cache_path)
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
        datasets = dsg.flatten().datasets
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
