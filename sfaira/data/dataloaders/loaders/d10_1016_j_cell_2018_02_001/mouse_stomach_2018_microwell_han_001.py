import os
from typing import Union
from .base import Dataset_d10_1016_j_cell_2018_02_001


class Dataset(Dataset_d10_1016_j_cell_2018_02_001):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "mouse_stomach_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.organ = "stomach"

        self.class_maps = {
            "0": {
                "Antral mucous cell (Stomach)": "antral mucous cell",
                "Dendritic cell(Stomach)": "dendritic cell",
                "Dividing cell(Stomach)": "proliferative cell",
                "Epithelial cell_Gkn3 high(Stomach)": "epithelial cell",
                "Epithelial cell_Krt20 high(Stomach)": "epithelial cell",
                "Epithelial cell_Pla2g1b high(Stomach)": "epithelial cell",
                "G cell(Stomach)": "G cell",
                "Gastric mucosal cell(Stomach)": "gastric mucosal cell",
                "Macrophage(Stomach)": "macrophage",
                "Muscle cell(Stomach)": "muscle cell",
                "Parietal cell (Stomach)": "parietal cell",
                "Pit cell_Gm26917 high(Stomach)": "pit cell",
                "Pit cell_Ifrd1 high(Stomach)": "pit cell",
                "Stomach cell_Gkn2 high(Stomach)": "stomach cell",
                "Stomach cell_Mt2 high(Stomach)": "stomach cell",
                "Stomach cell_Muc5ac high(Stomach)": "stomach cell",
                "Tuft cell(Stomach)": "tuft cell"
            },
        }

    def _load(self):
        self._load_generalized(samplename="Stomach_dge")
