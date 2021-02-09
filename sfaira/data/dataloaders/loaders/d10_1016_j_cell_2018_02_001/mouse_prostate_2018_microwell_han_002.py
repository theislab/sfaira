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
        self.id = "mouse_prostate_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.organ = "prostate"

        self.class_maps = {
            "0": {
                "Dendritic cell(Prostate)": "dendritic cell",
                "Epithelial cell(Prostate)": "epithelial cell",
                "Glandular epithelium(Prostate)": "glandular epithelial cell",
                "Prostate gland cell(Prostate)": "glandular cell",
                "Stromal cell(Prostate)": "stromal cell",
                "T cell(Prostate)": "T cell",
            },
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, samplename="Prostate2_dge")
