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
        self.id = "human_prostate_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = "prostate gland"
        self.class_maps = {
            "0": {
                "Antigen presenting cell (RPS high)": "Antigen presenting cell (RPS high)",
                "Basal cell": "Basal cell",
                "Dendritic cell": "Dendritic cell",
                "Endothelial cell": "Endothelial cell",
                "Endothelial cell (APC)": "Endothelial cell",
                "Endothelial cell (endothelial to mesenchymal transition)": "Endothelial cell",
                "Enterocyte progenitor": "Enterocyte progenitor",
                "Epithelial cell (intermediated)": "Epithelial cell (intermediated)",
                "Fasciculata cell": "Fasciculata cell",
                "Fetal enterocyte": "Fetal enterocyte",
                "Fetal epithelial progenitor": "Fetal epithelial progenitor",
                "Gastric endocrine cell": "Gastric endocrine cell",
                "Goblet cell": "Goblet cell",
                "Macrophage": "Macrophage",
                "Monocyte": "Monocyte",
                "Primordial germ cell": "Primordial germ cell",
                "Smooth muscle cell": "Smooth muscle cell",
                "Stratified epithelial cell": "Stratified epithelial cell",
                "Stromal cell": "Stromal cell",
                "T cell": "T cell",
            },
        }

    def _load(self):
        self._load_generalized(sample_id="AdultProstate_1")
