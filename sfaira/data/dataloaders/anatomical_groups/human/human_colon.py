from typing import Union

from sfaira.data.base import DatasetGroup
from sfaira.data.dataloaders.super_group import DatasetSuperGroupSfaira


class DatasetGroupColon(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(path=path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "human_colon_2019_10x_kinchen_001",
            "human_colon_2019_10x_smilie_001",
            "human_colon_2019_10x_wang_001",
            "human_colon_2020_10x_james_001",
            "human_colon_2020_microwell_han_001_10.1038/s41586-020-2157-4",
            "human_colon_2020_microwell_han_002_10.1038/s41586-020-2157-4",
            "human_colon_2020_microwell_han_003_10.1038/s41586-020-2157-4",
            "human_colon_2020_microwell_han_004_10.1038/s41586-020-2157-4",
        ])
        datasets = dsg.flatten().datasets
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
