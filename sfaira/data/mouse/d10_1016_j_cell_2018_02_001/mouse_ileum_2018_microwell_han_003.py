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
        self.id = "mouse_ileum_2018_microwell-seq_han_003_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "ileum"
        self.sub_tissue = "ileum"

        self.class_maps = {
            "0": {
                'B cell_Ighd high(Small-Intestine)': 'B cell',
                'B cell_Igkv12-46 high(Small-Intestine)': 'B cell',
                'B cell_Jchain high(Small-Intestine)': 'B cell',
                'B cell_Ms4a1 high(Small-Intestine)': 'B cell',
                'Columnar epithelium(Small-Intestine)': 'epithelial cell',
                'Dendritic cell_Siglech high(Small-Intestine)': 'dendritic cell',
                'Dendrtic cell_Cst3 high(Small-Intestine)': 'dendritic cell',
                'Epithelial cell_Kcne3 high(Small-Intestine)': 'epithelial cell',
                'Epithelial cell_Sh2d6 high(Small-Intestine)': 'epithelial cell',
                'Epithelium of small intestinal villi_Fabp1 high(Small-Intestine)': 'epithelial cell villi',
                'Epithelium of small intestinal villi_Fabp6 high(Small-Intestine)': 'epithelial cell villi',
                'Epithelium of small intestinal villi_Gm23935 high(Small-Intestine)': 'epithelial cell villi',
                'Epithelium of small intestinal villi_mt-Nd1 high(Small-Intestine)': 'epithelial cell villi',
                'Macrophage_Apoe high(Small-Intestine)': 'macrophage',
                'Macrophage_Cxcl2 high(Small-Intestine)': 'macrophage',
                'Paneth cell(Small-Intestine)': 'paneth cell',
                'S cell_Chgb high(Small-Intestine)': 'enteroendocrine cell',
                'S cell_Gip high(Small-Intestine)': 'enteroendocrine cell',
                'Stromal cell_Adamdec1 high(Small-Intestine)': 'stromal cell',
                'Stromal cell_Dcn high(Small-Intestine)': 'stromal cell',
                'T cell_Ccl5 high(Small-Intestine)': 'T cell',
                'T cell_Icos high(Small-Intestine)': 'T cell',
                'T cell_Cd7 high(Small-Intestine)': 'T cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "SmallIntestine3_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
