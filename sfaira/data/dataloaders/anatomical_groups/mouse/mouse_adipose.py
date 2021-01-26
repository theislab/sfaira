from typing import Union

from sfaira.data import DatasetGroup, DatasetSuperGroupSfaira


class DatasetGroupAdipose(DatasetGroup):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        dsg = DatasetSuperGroupSfaira(path=path, meta_path=meta_path, cache_path=cache_path)
        dsg.subset(key="id", values=[
            "mouse_adipose_2019_10x_pisco_001_10.1101/661728",
            "mouse_adipose_2019_smartseq2_pisco_001_10.1101/661728",
            "mouse_adipose_2019_smartseq2_pisco_002_10.1101/661728",
            "mouse_adipose_2019_smartseq2_pisco_003_10.1101/661728",
            "mouse_adipose_2019_smartseq2_pisco_004_10.1101/661728",
        ])
        datasets = dsg.flatten().datasets
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))
