from typing import Union
from .base import Dataset_d10_1038_s41586_020_2157_4


class Dataset(Dataset_d10_1038_s41586_020_2157_4):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_lung_2020_microwell_han_002_10.1038/s41586-020-2157-4"
        self.organ = "lung"
        self.class_maps = {
            "0": {
                "AT2 cell": "AT2",
                "Antigen presenting cell (RPS high)": "unknown",
                "B cell": "B cell lineage",
                "B cell (Plasmocyte)": "B cell lineage",
                "Basal cell": "Basal",
                "CB CD34+": "Fetal airway progenitors",
                "Chondrocyte": "1_Stroma",
                "Dendritic cell": "Dendritic cells",
                "Endothelial cell": "1_Endothelial",
                "Endothelial cell (APC)": "1_Endothelial",
                "Endothelial cell (endothelial to mesenchymal transition)": "1_Endothelial",
                "Enterocyte progenitor": "1_Epithelial",
                "Epithelial cell": "1_Epithelial",
                "Epithelial cell (intermediated)": "1_Epithelial",
                "Erythroid cell": "Erythrocytes",
                "Erythroid progenitor cell (RP high)": "Erythrocytes",
                "Fasciculata cell": "unknown",
                "Fetal Neuron": "unknown",
                "Fetal chondrocyte": "1_Stroma",
                "Fetal endocrine cell": "unknown",
                "Fetal enterocyte ": "1_Epithelial",
                "Fetal epithelial progenitor": "1_Epithelial",
                "Fetal fibroblast": "Fibroblasts",
                "Fetal mesenchymal progenitor": "1_Stroma",
                "Fetal neuron": "unknown",
                "Fetal skeletal muscle cell": "unknown",
                "Fetal stromal cell": "1_Stroma",
                "Fibroblast": "Fibroblasts",
                "Gastric endocrine cell": "unknown",
                "Goblet cell": "Secretory",
                "Kidney intercalated cell": "unknown",
                "Loop of Henle": "unknown",
                "M2 Macrophage": "Macrophages",
                "Macrophage": "Macrophages",
                "Mast cell": "Mast cells",
                "Mesothelial cell": "Mast cells",
                "Monocyte": "Monocytes",
                "Myeloid cell": "2_Myeloid",
                "Neutrophil": "Neutrophilic",
                "Neutrophil (RPS high)": "Neutrophilic",
                "Primordial germ cell": "unknown",
                "Proliferating T cell": "T cell lineage",
                "Proximal tubule progenitor": "unknown",
                "Sinusoidal endothelial cell": "1_Endothelial",
                "Smooth muscle cell": "2_Smooth Muscle",
                "Stratified epithelial cell": "1_Epithelial",
                "Stromal cell": "1_Stroma",
                "T cell": "T cell lineage",
                "Ventricle cardiomyocyte": "1_Stroma",
                "hESC": "Fetal airway progenitors",
            },
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, sample_id="AdultLung_3")
