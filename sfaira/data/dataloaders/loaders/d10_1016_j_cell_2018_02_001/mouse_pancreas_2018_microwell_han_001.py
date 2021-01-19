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
        self.id = "mouse_pancreas_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download_url_data = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "pancreas"

        self.class_maps = {
            "0": {
                "Acinar cell(Pancreas)": "pancreatic acinar cell",
                "Dendrtic cell(Pancreas)": "dendritic cell",
                "Ductal cell(Pancreas)": "pancreatic ductal cell",
                "Endocrine cell(Pancreas)": "endocrine cell",
                "Dividing cell(Pancreas)": "endocrine cell",
                "Endothelial cell_Fabp4 high(Pancreas)": "endothelial cell",
                "Endothelial cell_Lrg1 high(Pancreas)": "endothelial cell",
                "Endothelial cell_Tm4sf1 high(Pancreas)": "endothelial cell",
                "Erythroblast_Hbb-bt high(Pancreas)": "erythroblast",
                "Erythroblast_Igkc high(Pancreas)": "erythroblast",
                "Granulocyte(Pancreas)": "granulocyte",
                "Macrophage_Ly6c2 high(Pancreas)": "macrophage",
                "Macrophage(Pancreas)": "macrophage",
                "Glial cell(Pancreas)": "glial cell",
                "Smooth muscle cell_Acta2 high(Pancreas)": "smooth muscle cell",
                "Smooth muscle cell_Rgs5 high(Pancreas)": "smooth muscle cell",
                "Stromal cell_Fn1 high(Pancreas)": "stromal cell",
                "Stromal cell_Mfap4 high(Pancreas)": "stromal cell",
                "Stromal cell_Smoc2 high(Pancreas)": "stromal cell",
                "T cell(Pancreas)": "t cell",
                "B cell(Pancreas)": "b cell",
                "β-cell(Pancreas)": "pancreatic B cell"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Pancreas_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
