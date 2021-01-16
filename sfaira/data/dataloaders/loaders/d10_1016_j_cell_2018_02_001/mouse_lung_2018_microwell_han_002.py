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
        self.id = "mouse_lung_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "lung"

        self.class_maps = {
            "0": {
                "AT1 Cell(Lung)": "alveolar epithelial cell type I",
                "AT2 Cell(Lung)": "alveolar epithelial cell type II",
                "Alveolar bipotent progenitor(Lung)": "alveolar bipotent progenitor",
                "Alveolar macrophage_Ear2 high(Lung)": "alveolar macrophage",
                "Alveolar macrophage_Pclaf high(Lung)": "alveolar macrophage",
                "B Cell(Lung)": "B cell",
                "Basophil(Lung)": "basophil",
                "Ciliated cell(Lung)": "ciliated cell",
                "Clara Cell(Lung)": "clara cell",
                "Conventional dendritic cell_Gngt2 high(Lung)": "dendritic cell",
                "Conventional dendritic cell_H2-M2 high(Lung)": "dendritic cell",
                "Conventional dendritic cell_Mgl2 high(Lung)": "dendritic cell",
                "Conventional dendritic cell_Tubb5 high(Lung)": "dendritic cell",
                "Dendritic cell_Naaa high(Lung)": "dendritic cell",
                "Dividing T cells(Lung)": "T cell",
                "Dividing cells(Lung)": "unknown",
                "Dividing dendritic cells(Lung)": "dendritic cell",
                "Endothelial cell_Kdr high(Lung)": "endothelial cell",
                "Endothelial cell_Tmem100 high(Lung)": "endothelial cell",
                "Endothelial cells_Vwf high(Lung)": "endothelial cell",
                "Eosinophil granulocyte(Lung)": "eosinophil",
                "Ig−producing B cell(Lung)": "B cell",
                "Interstitial macrophage(Lung)": "lung macrophage",
                "Monocyte progenitor cell(Lung)": "monocyte progenitor",
                "NK Cell(Lung)": "NK cell",
                "Neutrophil granulocyte(Lung)": "neutrophil",
                "Nuocyte(Lung)": "nuocyte",
                "Plasmacytoid dendritic cell(Lung)": "plasmacytoid dendritic cell",
                "Stromal cell_Acta2 high(Lung)": "stromal cell",
                "Stromal cell_Dcn high(Lung)": "stromal cell",
                "Stromal cell_Inmt high(Lung)": "stromal cell",
                "T Cell_Cd8b1 high(Lung)": "CD8-positive, alpha-beta T cell",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Lung2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
