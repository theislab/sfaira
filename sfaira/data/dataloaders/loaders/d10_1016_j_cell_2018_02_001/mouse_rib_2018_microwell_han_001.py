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
        self.id = "mouse_rib_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.organ = "rib"

        self.class_maps = {
            "0": {
                "B cell(Neonatal-Rib)": "B cell",
                "Cartilage cell_Clu high(Neonatal-Rib)": "cartilage cell",
                "Cartilage cell_Col2a1 high(Neonatal-Rib)": "cartilage cell",
                "Cartilage cell_Cxcl14 high(Neonatal-Rib)": "cartilage cell",
                "Cartilage cell_Ppa1 high(Neonatal-Rib)": "cartilage cell",
                "Cartilage cell_Prg4 high(Neonatal-Rib)": "cartilage cell",
                "Dividing cell(Neonatal-Rib)": "proliferative cell",
                "Endothelial cell(Neonatal-Rib)": "endothelial cell",
                "Erythroblast_Hba-a1 high(Neonatal-Rib)": "erythroblast",
                "Erythroblast_Ttr high(Neonatal-Rib)": "erythroblast",
                "Granulocyte(Neonatal-Rib)": "granulocyte",
                "Macrophage_C1qc high(Neonatal-Rib)": "macrophage",
                "Macrophage_Ctss high(Neonatal-Rib)": "macrophage",
                "Muscle cell(Neonatal-Rib)": "muscle cell",
                "Muscle cell_Acta2 high(Neonatal-Rib)": "muscle cell",
                "Muscle cell_Actc1 high(Neonatal-Rib)": "muscle cell",
                "Neuron_Mpz high(Neonatal-Rib)": "neuron",
                "Neuron_Stmn2 high(Neonatal-Rib)": "neuron",
                "Neutrophil(Neonatal-Rib)": "neutrophil",
                "Neutrophil_Elane high(Neonatal-Rib)": "neutrophil",
                "Oligodendrocyte(Neonatal-Rib)": "oligodendrocyte",
                "Osteoblast(Neonatal-Rib)": "osteoblast",
                "Osteoclast(Neonatal-Rib)": "osteoclast",
                "Stromal cell_Acta1 high(Neonatal-Rib)": "stromal cell",
                "Stromal cell_Tnmd high(Neonatal-Rib)": "stromal cell",
            },
        }

    def _load(self):
        self._load_generalized(samplename="NeonatalRib1_dge")
