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
        self.id = "human_esophagus_2020_microwell_han_002_10.1038/s41586-020-2157-4"
        self.organ = "Esophagus"
        self.class_maps = {
            "0": {
                "Fibroblast": "Fibroblast",
                "Basal cell": "Basal cell",
                "Stratified epithelial cell": "Stratified epithelial cell",
                "Endothelial cell (APC)": "Endothelial cell (APC)",
                "Macrophage": "Macrophage",
                "B cell": "B cell",
                "T cell": "T cell",
                "Dendritic cell": "Dendritic cell",
                "Mast cell": "Mast cell",
                "B cell (Plasmocyte)": "B cell (Plasmocyte)",
                "Stromal cell": "Stromal cell",
                "Monocyte": "Monocyte",
                "Smooth muscle cell": "Smooth muscle cell",
                "Endothelial cell": "Endothelial cell",
                "Neutrophil": "Neutrophil",
                "Endothelial cell (endothelial to mesenchymal transition)": "Endothelial cell (endothelial to mesenchymal transition)",
                "Fetal stromal cell": "Fetal stromal cell",
                "CB CD34+": "CB CD34+",
                "Fetal epithelial progenitor": "Fetal epithelial progenitor",
                "Gastric endocrine cell": "Gastric endocrine cell",
                "Erythroid progenitor cell (RP high)": "Erythroid progenitor cell (RP high)",
                "Sinusoidal endothelial cell": "Sinusoidal endothelial cell",
                "Loop of Henle": "Loop of Henle",
                "Fetal mesenchymal progenitor": "Fetal mesenchymal progenitor",
            },
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, sample_id="AdultEsophagus_2")
