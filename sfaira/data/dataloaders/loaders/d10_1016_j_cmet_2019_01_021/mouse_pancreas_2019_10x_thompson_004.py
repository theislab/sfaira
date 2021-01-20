import os
from typing import Union
from .base import Dataset_d10_1016_j_cmet_2019_01_021


class Dataset(Dataset_d10_1016_j_cmet_2019_01_021):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "mouse_pancreas_2019_10x_thompson_004_10.1016/j.cmet.2019.01.021"

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "pancreas", "GSM3308549_NOD_14w_B")
            fn_meta = os.path.join(self.path, "mouse", "pancreas", "GSM3308549_NOD_14w_B_annotation.csv")
        else:
            fn_meta = os.path.join(fn, "_annotation.csv")

        data_path, sample_name = os.path.split(fn)
        self._load_generalized(data_path=data_path, sample_name=sample_name, path_meta=fn_meta)
