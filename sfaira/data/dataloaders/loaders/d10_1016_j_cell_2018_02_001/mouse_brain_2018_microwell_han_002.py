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
        self.id = "mouse_brain_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.organ = "brain"

        self.class_maps = {
            "0": {
                "Astroglial cell(Bergman glia)(Brain)": "Bergmann glial cell",
                "Astrocyte_Atp1b2 high(Brain)": "astrocyte",
                "Astrocyte_Mfe8 high(Brain)": "astrocyte",
                "Astrocyte_Pla2g7 high(Brain)": "astrocyte",
                "Granulocyte_Ngp high(Brain)": "granulocyte",
                "Hypothalamic ependymal cell(Brain)": "ependymal cell",
                "Macrophage_Klf2 high(Brain)": "macrophage",
                "Macrophage_Lyz2 high(Brain)": "macrophage",
                "Microglia(Brain)": "microglial cell",
                "Myelinating oligodendrocyte(Brain)": "oligodendrocyte",
                "Oligodendrocyte precursor cell(Brain)": "oligodendrocyte precursor cell",
                "Neuron(Brain)": "neuron",
                "Pan-GABAergic(Brain)": "GABAergic cell",
                "Schwann cell(Brain)": "schwann cell"
            },
        }

    def _load(self):
        self._load_generalized(samplename="Brain2_dge")
