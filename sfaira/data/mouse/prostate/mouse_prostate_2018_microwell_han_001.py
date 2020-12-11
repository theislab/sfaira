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

        self.id = "mouse_prostate_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "prostate"
        self.sub_tissue = "prostate"

        self.class_maps = {
            "0": {
                'Dendritic cell(Prostate)': 'dendritic cell',
                'Epithelial cell(Prostate)': 'epithelial cell',
                'Glandular epithelium(Prostate)': 'glandular epithelial cell',
                'Prostate gland cell(Prostate)': 'glandular cell',
                'Stromal cell(Prostate)': 'stromal cell',
                'T cell(Prostate)': 'T cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Prostate1_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
