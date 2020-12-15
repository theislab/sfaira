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
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.id = "mouse_mammarygland_2018_microwell-seq_han_004_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "mammarygland"
        self.sub_tissue = "mammarygland"

        self.class_maps = {
            "0": {
                'B cell_Cd79a&Fcer2a high(Mammary-Gland-Virgin)': 'B cell',
                'B cell_Cd79a&Iglc2 high(Mammary-Gland-Virgin)': 'B cell',
                'B cell_Jchain high(Mammary-Gland-Virgin)': 'B cell',
                'Dendritic cell_Cst3 high(Mammary-Gland-Virgin)': 'dendritic cell',
                'Dendritic cell_Fscn1 high(Mammary-Gland-Virgin)': 'dendritic cell',
                'Dendritic cell_Siglech high(Mammary-Gland-Virgin)': 'dendritic cell',
                'Dividing cell(Mammary-Gland-Virgin)': 'proliferative cell',
                'Luminal cell_Krt19 high (Mammary-Gland-Virgin)': 'luminal epithelial cell of mammary gland',
                'Luminal progenitor(Mammary-Gland-Virgin)': 'luminal progenitor cell',
                'Macrophage_C1qc high(Mammary-Gland-Virgin)': 'macrophage',
                'Macrophage_Lyz1 high(Mammary-Gland-Virgin)': 'macrophage',
                'NK cell(Mammary-Gland-Virgin)': 'NK cell',
                'Stem and progenitor cell(Mammary-Gland-Virgin)': 'stem and progenitor cell',
                'Stromal cell_Col3a1 high(Mammary-Gland-Virgin)': 'stromal cell',
                'Stromal cell_Pi16 high(Mammary-Gland-Virgin)': 'stromal cell',
                'T cell_Cd8b1 high(Mammary-Gland-Virgin)': 'T cell',
                'T cell_Ly6c2 high(Mammary-Gland-Virgin)': 'T cell',
                'T-cells_Ctla4 high(Mammary-Gland-Virgin)': 'T cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "MammaryGland.Virgin4_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
