from typing import Union

from sfaira.data.base import DatasetGroup
from sfaira.data.dataloaders.super_group import DatasetSuperGroupSfaira


class DatasetGroupKidney(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(path=path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "human_kidney_2019_10xSn_lake_001",
            "human_kidney_2019_10x_stewart_001",
            "human_kidney_2020_10x_liao_001",
            "human_kidney_2020_microwell_han_001_10.1038/s41586-020-2157-4",
            "human_kidney_2020_microwell_han_002_10.1038/s41586-020-2157-4",
            "human_kidney_2020_microwell_han_003_10.1038/s41586-020-2157-4",
            "human_kidney_2020_microwell_han_004_10.1038/s41586-020-2157-4",
            "human_kidney_2020_microwell_han_005_10.1038/s41586-020-2157-4",
            "human_kidney_2020_microwell_han_006_10.1038/s41586-020-2157-4",
            "human_kidney_2020_microwell_han_007_10.1038/s41586-020-2157-4",
        ])
        super().__init__(datasets=dsg.flatten().datasets)
