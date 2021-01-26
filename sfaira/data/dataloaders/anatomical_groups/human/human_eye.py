from typing import Union

from sfaira.data import DatasetGroup, DatasetSuperGroupSfaira


class DatasetGroupEye(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(path=path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "human_eye_2019_10x_lukowski_001",
            "human_eye_2019_10x_menon_001",
            "human_eye_2019_10x_voigt_001",
            "human_eye_2020_microwell_han_001_10.1038/s41586-020-2157-4",
        ])
        datasets = dsg.flatten().datasets
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
