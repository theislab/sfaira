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
        self.id = "mouse_uterus_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "uterus"

        self.class_maps = {
            "0": {
                'B cell(Uterus)': 'B cell',
                'Dendritic cell(Uterus)': 'dendritic cell',
                'Endothelial cell_Cldn5 high(Uterus)': 'endothelial cell',
                'Endothelial cell_Tm4sf1 high(Uterus)': 'endothelial cell',
                'Glandular epithelium_Ltf high(Uterus)': 'glandular epithelial cell',
                'Glandular epithelium_Sprr2f high(Uterus)': 'glandular epithelial cell',
                'Granulocyte(Uterus)': 'granulocyte',
                'Keratinocyte(Uterus)': 'keratinocyte',
                'Macrophage(Uterus)': 'macrophage',
                'Monocyte(Uterus)': 'monocyte',
                'Muscle cell_Mgp high(Uterus)': 'muscle cell',
                'Muscle cell_Pcp4 high(Uterus)': 'muscle cell',
                'Smooth muscle cell_Rgs5 high(Uterus)': 'smooth muscle cell',
                'NK cell(Uterus)': 'NK cell',
                'Stromal cell_Ccl11 high(Uterus)': 'stromal cell',
                'Stromal cell_Cxcl14 high(Uterus)': 'stromal cell',
                'Stromal cell_Gm23935 high(Uterus)': 'stromal cell',
                'Stromal cell_Has1 high(Uterus)': 'stromal cell',
                'Stromal cell_Hsd11b2 high(Uterus)': 'stromal cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Uterus1_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
