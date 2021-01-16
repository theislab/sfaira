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
        self.id = "mouse_femalegonad_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "femalegonad"

        self.class_maps = {
            "0": {
                "Cumulus cell_Car14 high(Ovary)": "cumulus cell",
                "Cumulus cell_Nupr1 high(Ovary)": "cumulus cell",
                "Cumulus cell_Ube2c high(Ovary)": "cumulus cell",
                "Granulosa cell_Inhba high(Ovary)": "granulosa cell",
                "Granulosa cell_Kctd14 high(Ovary)": "granulosa cell",
                "Large luteal cell(Ovary)": "large luteal cell",
                "Macrophage_Lyz2 high(Ovary)": "macrophage",
                "Marcrophage_Cd74 high(Ovary)": "macrophage",
                "Ovarian surface epithelium cell(Ovary)": "epithelial cell of ovarian surface",
                "Ovarian vascular surface endothelium cell(Ovary)": "endothelial cell of ovarian surface",
                "Small luteal cell(Ovary)": "small luteal cell",
                "Stroma cell (Ovary)": "stromal cell",
                "Thecal cell(Ovary)": "thecal cell",
                "luteal cells(Ovary)": "luteal cell",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Ovary1_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
