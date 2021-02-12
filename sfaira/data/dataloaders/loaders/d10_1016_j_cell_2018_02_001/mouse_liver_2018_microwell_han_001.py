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
        self.id = "mouse_liver_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.organ = "liver"

        self.class_maps = {
            "0": {
                "B cell_Fcmr high(Liver)": "B cell",
                "B cell_Jchain high(Liver)": "B cell",
                "Dendritic cell_Cst3 high(Liver)": "dendritic cell",
                "Dendritic cell_Siglech high(Liver)": "dendritic cell",
                "Endothelial cell(Liver)": "endothelial cell of hepatic sinusoid",
                "Epithelial cell(Liver)": "duct epithelial cell",
                "Epithelia cell_Spp1 high(Liver)": "duct epithelial cell",
                "Erythroblast_Hbb-bs high(Liver)": "erythroblast",
                "Erythroblast_Hbb-bt high(Liver)": "erythroblast",
                "Granulocyte(Liver)": "granulocyte",
                "Hepatocyte_Fabp1 high(Liver)": "hepatocyte",
                "Hepatocyte_mt-Nd4 high(Liver)": "hepatocyte",
                "Pericentral (PC) hepatocytes(Liver)": "hepatocyte",
                "Periportal (PP) hepatocyte(Liver)": "hepatocyte",
                "Kuppfer cell(Liver)": "Kupffer cell",
                "Macrophage_Chil3 high(Liver)": "macrophage",
                "Neutrophil_Ngp high(Liver)": "neutrophil",
                "Stromal cell(Liver)": "stromal cell",
                "T cell_Gzma high(Liver)": "T cell",
                "T cell_Trbc2 high(Liver)": "T cell",
            },
        }

    def _load(self):
        self._load_generalized(samplename="Liver1_dge")
