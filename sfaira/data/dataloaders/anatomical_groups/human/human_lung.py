from typing import Union

from sfaira.data.base import DatasetGroup
from sfaira.data.dataloaders.super_group import DatasetSuperGroupSfaira


class DatasetGroupLung(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(path=path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "human_lung_2019_10x_braga_001",
            "human_lung_2019_10x_braga_002",
            "human_lung_2019_dropseq_braga_003",
            "human_lung_2019_10x_madissoon_001",
            "human_lung_2020_10x_habermann_001",
            "human_lung_2020_10x_lukassen_001",
            "human_lung_2020_10x_lukassen_002",
            "human_lung_2020_10x_miller_001",
            "human_lung_2020_10x_travaglini_001",
            "human_lung_2020_microwell_han_001_10.1038/s41586-020-2157-4",
            "human_lung_2020_microwell_han_002_10.1038/s41586-020-2157-4",
            "human_lung_2020_microwell_han_003_10.1038/s41586-020-2157-4",
            "human_lung_2020_microwell_han_004_10.1038/s41586-020-2157-4",
            "human_lung_2020_microwell_han_005_10.1038/s41586-020-2157-4",
            "human_lung_2020_smartseq2_travaglini_002",
        ])
        super().__init__(datasets=dsg.flatten().datasets)
