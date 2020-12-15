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
        self.id = "mouse_spleen_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "spleen"
        self.sub_tissue = "spleen"

        self.class_maps = {
            "0": {
                'Erythroblast(Spleen)': 'proerythroblast',
                'Dendritic cell_S100a4 high(Spleen)': 'dendritic cell',
                'Dendritic cell_Siglech high(Spleen)': 'dendritic cell',
                'Granulocyte(Spleen)': 'granulocyte',
                'Macrophage(Spleen)': 'macrophage',
                'Monocyte(Spleen)': 'monocyte',
                'NK cell(Spleen)': 'NK cell',
                'Neutrophil(Spleen)': 'neutrophil',
                'Plasma cell(Spleen)': 'plasma cell',
                'T cell(Spleen)': 'T cell',
                'Marginal zone B cell(Spleen)': 'B cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Spleen_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
        self.set_unkown_class_id(ids=[np.nan, "nan"])

        