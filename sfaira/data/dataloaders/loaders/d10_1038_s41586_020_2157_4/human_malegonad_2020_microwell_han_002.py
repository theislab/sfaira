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
        self.id = "human_testis_2020_microwell_han_002_10.1038/s41586-020-2157-4"
        self.organ = "testis"
        self.class_maps = {
            "0": {
                "Antigen presenting cell (RPS high)": "Antigen presenting cell (RPS high)",
                "B cell": "B cell",
                "CB CD34+": "CB CD34+",
                "Dendritic cell": "Dendritic cell",
                "Endothelial cell": "Endothelial cells",
                "Erythroid cell": "Erythroid cell",
                "Erythroid progenitor cell (RP high)": "Erythroid progenitor cell (RP high)",
                "Fasciculata cell": "Fasciculata cell",
                "Fetal acinar cell": "Fetal acinar cell",
                "Fetal chondrocyte": "Fetal chondrocyte",
                "Fetal epithelial progenitor": "Fetal epithelial progenitor",
                "Fetal fibroblast": "Fetal fibroblast",
                "Fetal mesenchymal progenitor": "Fetal mesenchymal progenitor",
                "Fetal neuron": "Fetal neuron",
                "Fetal skeletal muscle cell": "Fetal skeletal muscle cell",
                "Fetal stromal cell": "Fetal stromal cell",
                "Immature sertoli cell (Pre-Sertoli cell)": "Sertoli cells",
                "Loop of Henle": "Loop of Henle",
                "Macrophage": "Macrophages",
                "Monocyte": "Monocyte",
                "Neutrophil": "Neutrophil",
                "Neutrophil (RPS high)": "Neutrophil (RPS high)",
                "Primordial germ cell": "Primordial germ cell",
                "Proximal tubule progenitor": "Proximal tubule progenitor",
                "Smooth muscle cell": "Smooth muscle cell",
                "Stromal cell": "Stromal cell",
                "T cell": "T cell",
                "Ureteric bud cell": "Ureteric bud cell",
            },
        }

    def _load(self):
        self._load_generalized(sample_id="FetalMaleGonad_2")
