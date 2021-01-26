from typing import Union

from sfaira.data.base import DatasetGroup
from sfaira.data.dataloaders.super_group import DatasetSuperGroupSfaira


class DatasetGroupIleum(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(path=path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "mouse_ileum_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001",
            "mouse_ileum_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001",
            "mouse_ileum_2018_microwell-seq_han_003_10.1016/j.cell.2018.02.001",
        ])
        datasets = dsg.flatten().datasets
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
