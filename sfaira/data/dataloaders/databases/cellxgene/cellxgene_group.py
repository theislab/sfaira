import os
from typing import Union

from sfaira.data import DatasetGroup
from sfaira.consts import ADATA_IDS_CELLXGENE

from .cellxgene_loader import Dataset


class DatasetGroup(DatasetGroup):

    def __init__(
        self,
        data_path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        fn_ls = os.listdir(data_path)
        fn_ls = [x for x in fn_ls if x in ADATA_IDS_CELLXGENE.accepted_file_names]
        datasets = [
            Dataset(data_path=data_path, fn=x, meta_path=meta_path, cache_path=cache_path)
            for x in fn_ls
        ]
        keys = [x.id for x in datasets]
        super().__init__(dict(zip(keys, datasets)))
