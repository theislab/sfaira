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
        self.id = "mouse_malegonad_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download_url_data = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "malegonad"

        self.class_maps = {
            "0": {
                "Elongating spermatid(Testis)": "elongating spermatid",
                "Erythroblast_Hbb-bs high(Testis)": "erythroblast",
                "Leydig cell(Testis)": "leydig cell",
                "Macrophage_Lyz2 high(Testis)": "macrophage",
                "Pre-Sertoli cell_Cst9 high(Testis)": "pre-sertoli cell",
                "Pre-Sertoli cell_Ctsl high(Testis)": "pre-sertoli cell",
                "Preleptotene spermatogonia(Testis)": "preleptotene spermatogonia",
                "Sertoli cell(Testis)": "sertoli cell",
                "Spermatids_1700016P04Rik high(Testis)": "spermatid",
                "Spermatids_Cst13 high(Testis)": "spermatid",
                "Spermatids_Hmgb4 high(Testis)": "spermatid",
                "Spermatids_Tnp1 high(Testis)": "spermatid",
                "Spermatocyte_1700001F09Rik high(Testis)": "spermatocyte",
                "Spermatocyte_Cabs1 high(Testis)": "spermatocyte",
                "Spermatocyte_Calm2 high(Testis)": "spermatocyte",
                "Spermatocyte_Mesp1 high(Testis)": "spermatocyte",
                "Spermatocyte_Slc2a3 high(Testis)": "spermatocyte",
                "Spermatogonia_1700001P01Rik high(Testis)": "spermatogonia",
                "Spermatogonia_Tbc1d23 high(Testis)": "spermatogonia"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Testis2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
