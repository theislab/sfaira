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

        self.id = "mouse_liver_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download_url_data = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "liver"

        self.class_maps = {
            "0": {
                "B cell_Jchain high(Liver)": "B cell",
                "Dendritic cell_Cst3 high(Liver)": "dendritic cell",
                "Dendritic cell_Siglech high(Liver)": "dendritic cell",
                "Epithelial cell(Liver)": "duct epithelial cell",
                "Epithelia cell_Spp1 high(Liver)": "duct epithelial cell",
                "Erythroblast_Hbb-bs high(Liver)": "erythroblast",
                "Hepatocyte_Fabp1 high(Liver)": "hepatocyte",
                "Pericentral (PC) hepatocytes(Liver)": "hepatocyte",
                "Periportal (PP) hepatocyte(Liver)": "hepatocyte",
                "Kuppfer cell(Liver)": "Kupffer cell",
                "Macrophage_Chil3 high(Liver)": "macrophage",
                "Stromal cell(Liver)": "stromal cell",
                "T cell_Gzma high(Liver)": "T cell",
                "T cell_Trbc2 high(Liver)": "T cell",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Liver2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
