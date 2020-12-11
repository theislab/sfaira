import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetMca


class Dataset(DatasetMca):

    id: str

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetMca.__init__(self=self, path=path, meta_path=meta_path, **kwargs)

        self.id = "mouse_bladder_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "bladder"
        self.sub_tissue = "bladder"

        self.class_maps = {
            "0": {
                "Endothelial cell_Ly6c1 high(Bladder)": 'endothelial cell',
                "Vascular endothelial cell(Bladder)": 'endothelial cell',
                'Urothelium(Bladder)': 'bladder urothelial cell',
                'Dendritic cell_Cd74 high(Bladder)': 'dendritic cell',
                'Dendritic cell_Lyz2 high(Bladder)': 'dendritic cell',
                'Macrophage_Pf4 high(Bladder)': 'macrophage',
                'NK cell(Bladder)': 'NK cell',
                'Basal epithelial cell(Bladder)': 'basal epithelial cell',
                'Epithelial cell_Upk3a high(Bladder)': 'epithelial cell',
                'Epithelial cell_Gm23935 high(Bladder)': 'epithelial cell',
                'Mesenchymal stromal cell(Bladder)': 'mesenchymal stromal cell',
                'Stromal cell_Dpt high(Bladder)': 'stromal cell',
                'Stromal cell_Car3 high(Bladder)': 'stromal cell',
                'Smooth muscle cell(Bladder)': 'smooth muscle cell',
                'Vascular smooth muscle progenitor cell(Bladder)': 'smooth muscle cell',
                'Umbrella cell(Bladder)': 'umbrella cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Bladder_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
