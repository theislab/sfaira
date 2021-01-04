import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase


class Dataset_d10_1016_j_cmet_2019_01_021(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.author = "Bhushan"
        self.doi = "10.1016/j.cmet.2019.01.021"
        self.download = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117770"
        self.healthy = False
        self.normalization = 'raw'
        self.organ = "pancreas"
        self.protocol = "10x"
        self.state_exact = "diabetic"
        self.sub_tissue = "pancreas"
        self.year = 2019

        self.var_symbol_col = "index"
