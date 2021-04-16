from typing import Union

from sfaira.data import DatasetGroup
from sfaira.data.dataloaders.super_group import Universe


class DatasetGroupBrain(DatasetGroup):

    def __init__(
        self,
        data_path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = Universe(data_path=data_path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "human_brain_2017_DroNcSeq_habib_001",
            "human_brain_2020_microwell_han_001_10.1038/s41586-020-2157-4",
            "human_brain_2020_microwell_han_002_10.1038/s41586-020-2157-4",
            "human_brain_2020_microwell_han_003_10.1038/s41586-020-2157-4",
            "human_brain_2020_microwell_han_004_10.1038/s41586-020-2157-4",
            "human_brain_2020_microwell_han_005_10.1038/s41586-020-2157-4",
            "human_brain_2020_microwell_han_006_10.1038/s41586-020-2157-4",
        ])
        super().__init__(datasets=dsg.flatten().datasets)
