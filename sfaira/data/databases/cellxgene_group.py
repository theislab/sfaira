import os
from typing import Union

from .external import DatasetGroupBase

from .cellxgene_loader import Dataset


class DatasetGroupCellxgene(DatasetGroupBase):

    def __init__(
        self,
        path: Union[str, None] = None,
        meta_path: Union[str, None] = None
    ):
        fn_ls = os.listdir(path)
        fn_ls = [x for x in fn_ls if x in self.accepted_file_names]
        datasets = [
            Dataset(path=path, fn=x, meta_path=meta_path)
            for x in fn_ls
        ]
        keys = [x.id for x in datasets]
        self.datasets = dict(zip(keys, datasets))

    @property
    def accepted_file_names(self):
        return [
            "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remixed.h5ad"
        ]
