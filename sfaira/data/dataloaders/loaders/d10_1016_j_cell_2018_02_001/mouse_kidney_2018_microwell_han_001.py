import os
from typing import Union
from .base import Dataset_d10_1016_j_cell_2018_02_001


class Dataset(Dataset_d10_1016_j_cell_2018_02_001):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "mouse_kidney_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.organ = "kidney"

        self.class_maps = {
            "0": {
                "Cell in cell cycle(Fetal_Kidney)": "fetal proliferative cell",
                "Metanephric mesenchyme(Fetal_Kidney)": "fetal mesenchymal cell"
            },
        }

    def _load(self):
        self._load_generalized(samplename="Kidney1_dge")
