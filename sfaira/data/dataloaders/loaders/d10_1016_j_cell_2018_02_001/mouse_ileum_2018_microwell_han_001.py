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
        self.id = "mouse_ileum_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.organ = "ileum"

        self.class_maps = {
            "0": {
                "B cell_Ighd high(Small-Intestine)": "B cell",
                "B cell_Igkv12-46 high(Small-Intestine)": "B cell",
                "B cell_Jchain high(Small-Intestine)": "B cell",
                "B cell_Ms4a1 high(Small-Intestine)": "B cell",
                "Columnar epithelium(Small-Intestine)": "epithelial cell",
                "Dendritic cell_Siglech high(Small-Intestine)": "dendritic cell",
                "Dendrtic cell_Cst3 high(Small-Intestine)": "dendritic cell",
                "Epithelial cell_Kcne3 high(Small-Intestine)": "epithelial cell",
                "Epithelial cell_Sh2d6 high(Small-Intestine)": "epithelial cell",
                "Epithelium of small intestinal villi_Fabp1 high(Small-Intestine)": "epithelial cell villi",
                "Epithelium of small intestinal villi_Fabp6 high(Small-Intestine)": "epithelial cell villi",
                "Epithelium of small intestinal villi_Gm23935 high(Small-Intestine)": "epithelial cell villi",
                "Epithelium of small intestinal villi_mt-Nd1 high(Small-Intestine)": "epithelial cell villi",
                "Macrophage_Apoe high(Small-Intestine)": "macrophage",
                "Macrophage_Cxcl2 high(Small-Intestine)": "macrophage",
                "Paneth cell(Small-Intestine)": "paneth cell",
                "S cell_Chgb high(Small-Intestine)": "enteroendocrine cell",
                "S cell_Gip high(Small-Intestine)": "enteroendocrine cell",
                "Stromal cell_Adamdec1 high(Small-Intestine)": "stromal cell",
                "Stromal cell_Dcn high(Small-Intestine)": "stromal cell",
                "T cell_Ccl5 high(Small-Intestine)": "T cell",
                "T cell_Icos high(Small-Intestine)": "T cell",
                "T cell_Cd7 high(Small-Intestine)": "T cell",
            },
        }

    def _load(self):
        self._load_generalized(samplename="SmallIntestine1_dge")
