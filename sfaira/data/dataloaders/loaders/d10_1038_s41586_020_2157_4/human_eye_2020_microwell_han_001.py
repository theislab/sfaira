from typing import Union
from .base import Dataset_d10_1038_s41586_020_2157_4


class Dataset(Dataset_d10_1038_s41586_020_2157_4):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_eye_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = "Eye"
        self.class_maps = {
            "0": {
                "Fetal neuron": "Fetal neuron",
                "Fetal mesenchymal progenitor": "Fetal mesenchymal progenitor",
                "Fetal epithelial progenitor": "Fetal epithelial progenitor",
                "Erythroid cell": "Erythroid cell",
                "Primordial germ cell": "Primordial germ cell",
                "Endothelial cell": "Endothelial cell",
                "Fetal skeletal muscle cell": "Fetal skeletal muscle cell",
                "Fetal stromal cell": "Fetal stromal cell",
                "Fetal fibroblast": "Fibroblast",
                "Fetal Neuron": "Fetal neuron",
                "Erythroid progenitor cell (RP high)": "Erythroid progenitor cell (RP high)",
                "Antigen presenting cell (RPS high)": "Antigen presenting cell (RPS high)",
                "Dendritic cell": "Dendritic cell",
                "Fetal endocrine cell": "Fetal endocrine cell",
                "Macrophage": "Macrophage",
                "T cell": "T cell",
                "Basal cell": "Basal cell",
                "Gastric endocrine cell": "Gastric endocrine cell",
                "Goblet cell": "Goblet cell",
                "Epithelial cell (intermediated)": "Epithelial cell (intermediated)",
                "Stratified epithelial cell": "Stratified epithelial cell",
                "CB CD34+": "CB CD34_pos",
                "hESC": "hESC"
            },
        }

    def _load(self):
        self._load_generalized(sample_id="FetalEyes_1")
