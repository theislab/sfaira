import os
from typing import Union

from sfaira.data import DatasetGroup
from sfaira.consts import AdataIdsCellxgene

from .cellxgene_loader import Dataset


class DatasetGroupCellxgene(DatasetGroup):

    def __init__(
        self,
        data_path: Union[str, None] = None,
        meta_path: Union[str, None] = None,
        cache_path: Union[str, None] = None
    ):
        self._adata_ids_cellxgene = AdataIdsCellxgene()
        fn_ls = os.listdir(data_path)
        fn_ls = [x for x in fn_ls if x in self._adata_ids_cellxgene.accepted_file_names]
        datasets = [
            Dataset(data_path=data_path, fn=x, meta_path=meta_path, cache_path=cache_path)
            for x in fn_ls
        ]
        keys = [x.id for x in datasets]
        super().__init__(dict(zip(keys, datasets)))
