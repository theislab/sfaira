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
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "mouse_muscle_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "muscle"
        self.sub_tissue = "muscle"

        self.class_maps = {
            "0": {
                'B cell_Jchain high(Muscle)': 'B cell',
                'B cell_Vpreb3 high(Muscle)': 'B cell',
                'Dendritic cell(Muscle)': 'dendritic cell',
                'Endothelial cell(Muscle)': 'endothelial cell',
                'Erythroblast_Car1 high(Muscle)': 'erythroblast',
                'Erythroblast_Car2 high(Muscle)': 'erythroblast',
                'Granulocyte monocyte progenitor cell(Muscle)': 'monocyte progenitor',
                'Macrophage_Ms4a6c high(Muscle)': 'macrophage',
                'Macrophage_Retnla high(Muscle)': 'macrophage',
                'Muscle cell_Tnnc1 high(Muscle)': 'muscle cell',
                'Muscle cell_Tnnc2 high(Muscle)': 'muscle cell',
                'Muscle progenitor cell(Muscle)': 'skeletal muscle satellite cell',
                'Neutrophil_Camp high(Muscle)': 'neutrophil',
                'Neutrophil_Prg2 high(Muscle)': 'neutrophil',
                'Neutrophil_Retnlg high(Muscle)': 'neutrophil',
                'Stromal cell(Muscle)': 'stromal cell',
                'T cell(Muscle)': 'T cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Muscle_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
