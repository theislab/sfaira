import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetMca


class Dataset(DatasetMca):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetMca.__init__(self=self, path=path, meta_path=meta_path, **kwargs)

        self.id = "mouse_stomach_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "stomach"
        self.sub_tissue = "stomach"

        self.class_maps = {
            "0": {
                'Antral mucous cell (Stomach)': 'antral mucous cell',
                'Dendritic cell(Stomach)': 'dendritic cell',
                'Dividing cell(Stomach)': 'proliferative cell',
                'Epithelial cell_Gkn3 high(Stomach)': 'epithelial cell',
                'Epithelial cell_Krt20 high(Stomach)': 'epithelial cell',
                'Epithelial cell_Pla2g1b high(Stomach)': 'epithelial cell',
                'G cell(Stomach)': 'G cell',
                'Gastric mucosal cell(Stomach)': 'gastric mucosal cell',
                'Macrophage(Stomach)': 'macrophage',
                'Muscle cell(Stomach)': 'muscle cell',
                'Parietal cell (Stomach)': 'parietal cell',
                'Pit cell_Gm26917 high(Stomach)': 'pit cell',
                'Pit cell_Ifrd1 high(Stomach)': 'pit cell',
                'Stomach cell_Gkn2 high(Stomach)': 'stomach cell',
                'Stomach cell_Mt2 high(Stomach)': 'stomach cell',
                'Stomach cell_Muc5ac high(Stomach)': 'stomach cell',
                'Tuft cell(Stomach)': 'tuft cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas", "500more_dge", "Stomach_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
