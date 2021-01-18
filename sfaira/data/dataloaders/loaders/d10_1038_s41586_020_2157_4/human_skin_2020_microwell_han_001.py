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
        self.id = "human_skin_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = "skin"
        self.class_maps = {
            "0": {
                "Antigen presenting cell (RPS high)": "Antigen presenting cell (RPS high)",
                "B cell": "B cell",
                "Basal cell": "Basal cell",
                "CB CD34+": "CB CD34+",
                "Dendritic cell": "Dendritic cell",
                "Endothelial cell": "Endothelial cell",
                "Endothelial cell (APC)": "Endothelial cell (APC)",
                "Epithelial cell": "Epithelial cell",
                "Erythroid cell": "Erythroid cell",
                "Erythroid progenitor cell (RP high)": "Erythroid progenitor cell (RP high)",
                "Fetal Neuron": "Fetal Neuron",
                "Fetal epithelial progenitor": "Fetal epithelial progenitor",
                "Fetal fibroblast": "Fetal fibroblast",
                "Fetal mesenchymal progenitor": "Fetal mesenchymal progenitor",
                "Fetal skeletal muscle cell": "Fetal skeletal muscle cell",
                "Fetal stromal cell": "Fetal stromal cell",
                "Fibroblast": "Fibroblast",
                "Kidney intercalated cell": "Kidney intercalated cell",
                "Macrophage": "Macrophage",
                "Mast cell": "Mast cell",
                "Monocyte": "Monocyte",
                "Neutrophil": "Neutrophil",
                "Neutrophil (RPS high)": "Neutrophil (RPS high)",
                "Primordial germ cell": "Primordial germ cell",
                "Proliferating T cell": "Proliferating T cell",
                "Smooth muscle cell": "Smooth muscle cell",
                "Stromal cell": "Stromal cell",
                "T cell": "T cell",
                "hESC": "hESC",
            },
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, sample_id="FetalSkin_2")
