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
        self.id = "mouse_placenta_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.organ = "placenta"

        self.class_maps = {
            "0": {
                "B cell(Placenta)": "B cell",
                "Basophil(Placenta)": "basophil",
                "Decidual stromal cell(Placenta)": "decidual stromal cell",
                "Dendritic cell(Placenta)": "dendritic cell",
                "Endodermal cell_Afp high(Placenta)": "endodermal cell",
                "Endothelial cell_Maged2 high(Placenta)": "endothelial cell",
                "Erythroblast_Hbb-y high(Placenta)": "erythroblast",
                "Granulocyte monocyte progenitors(Placenta)": "monocyte progenitor",
                "Granulocyte_Neat1 high(Placenta)": "granulocyte",
                "Granulocyte_S100a9 high(Placenta)": "granulocyte",
                "HSPC_Lmo2 high(Placenta)": "HSPC",
                "Invasive spongiotrophoblast(Placenta)": "invasive spongiotrophoblast",
                "Labyrinthine trophoblast(Placenta)": "labyrinthine trophoblast",
                "Macrophage_Apoe high(Placenta)": "macrophage",
                "Macrophage_Spp1 high(Placenta)": "macrophage",
                "Megakaryocyte progenitor cell(Placenta)": "megakaryocte",
                "Monocyte(Placenta)": "monocyte",
                "NK cell(Placenta)": "NK cell",
                "NKT cell(Placenta)": "NKT cell",
                "PE lineage cell_Gkn2 high(Placenta)": "PE lineage cell",
                "PE lineage cell_S100g high(Placenta)": "PE lineage cell",
                "Progenitor trophoblast_Gjb3 high(Placenta)": "trophoblast progenitor",
                "Spiral artery trophoblast giant cells(Placenta)": "spiral artery trophoblast giant cells",
                "Spongiotrophoblast_Hsd11b2 high(Placenta)": "spongiotrophoblast",
                "Spongiotrophoblast_Phlda2 high(Placenta)": "spongiotrophoblast",
                "Stromal cell(Placenta)": "stromal cell",
                "Stromal cell_Acta2 high(Placenta)": "stromal cell",
                "Trophoblast progenitor_Taf7l high(Placenta)": "trophoblast progenitor",
            },
        }

    def _load(self):
        self._load_generalized(samplename="PlacentaE14.1_dge")
