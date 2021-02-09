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
        self.id = "mouse_spleen_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.organ = "spleen"

        self.class_maps = {
            "0": {
                "Erythroblast(Spleen)": "proerythroblast",
                "Dendritic cell_S100a4 high(Spleen)": "dendritic cell",
                "Dendritic cell_Siglech high(Spleen)": "dendritic cell",
                "Granulocyte(Spleen)": "granulocyte",
                "Macrophage(Spleen)": "macrophage",
                "Monocyte(Spleen)": "monocyte",
                "NK cell(Spleen)": "NK cell",
                "Neutrophil(Spleen)": "neutrophil",
                "Plasma cell(Spleen)": "plasma cell",
                "T cell(Spleen)": "T cell",
                "Marginal zone B cell(Spleen)": "B cell"
            },
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, samplename="Spleen_dge")
