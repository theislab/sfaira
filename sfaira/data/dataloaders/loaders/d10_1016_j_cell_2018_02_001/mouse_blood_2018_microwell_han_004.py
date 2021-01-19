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
        self.id = "mouse_blood_2018_microwell-seq_han_004_10.1016/j.cell.2018.02.001"
        self.download_url_data = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "blood"

        self.class_maps = {
            "0": {
                "B cell_Igha high(Peripheral_Blood)": "B cell",
                "B cell_Ly6d high(Peripheral_Blood)": "B cell",
                "B cell_Rps27rt high(Peripheral_Blood)": "B cell",
                "B cell_Vpreb3 high(Peripheral_Blood)": "B cell",
                "Basophil_Prss34 high(Peripheral_Blood)": "basophil",
                "Dendritic cell_Siglech high(Peripheral_Blood)": "dendritic cell",
                "Erythroblast_Car2 high(Peripheral_Blood)": "erythroblast",
                "Erythroblast_Hba-a2 high(Peripheral_Blood)": "erythroblast",
                "Macrophage_Ace high(Peripheral_Blood)": "macrophage",
                "Macrophage_Flt-ps1 high(Peripheral_Blood)": "macrophage",
                "Macrophage_Pf4 high(Peripheral_Blood)": "macrophage",
                "Macrophage_S100a4 high(Peripheral_Blood)": "macrophage",
                "Monocyte_Elane high(Peripheral_Blood)": "monocyte",
                "Monocyte_F13a1 high(Peripheral_Blood)": "monocyte",
                "NK cell_Gzma high(Peripheral_Blood)": "NK cell",
                "Neutrophil_Camp high(Peripheral_Blood)": "neutrophil",
                "Neutrophil_Il1b high(Peripheral_Blood)": "neutrophil",
                "Neutrophil_Ltf high(Peripheral_Blood)": "neutrophil",
                "Neutrophil_Retnlg high(Peripheral_Blood)": "neutrophil",
                "T cell_Gm14303 high(Peripheral_Blood)": "T cell",
                "T cell_Trbc2 high(Peripheral_Blood)": "T cell"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "PeripheralBlood4_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
