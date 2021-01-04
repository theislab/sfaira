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
        self.id = "mouse_thymus_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "thymus"
        self.sub_tissue = "thymus"

        self.class_maps = {
            "0": {
                'abT cell(Thymus)': 'abT cell',
                'B cell(Thymus)': "B cell",
                'DPT cell(Thymus)': "double positive T cell",
                'gdT cell (Thymus)': 'gdT cell',
                'Pre T cell(Thymus)': 'immature T cell',
                'Proliferating thymocyte(Thymus)': "immature T cell",
                'T cell_Id2 high(Thymus)': 'abT cell',  # TODO check, not sure about this gene
                'T cell_Ms4a4b high(Thymus)': 'abT cell'  # TODO check, not sure about this gene
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Thymus1_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
