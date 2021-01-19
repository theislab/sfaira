import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_ileum_2019_10x_wang_001_10.1084/jem.20191130"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/wang20_ileum.processed.h5ad"
        self.download_url_meta = None

        self.author = "Chen"
        self.doi = "10.1084/jem.20191130"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "ileum"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "Progenitor": "Progenitors",
                "Goblet": "Goblet cells",
                "Enterocyte": "Enterocytes",
                "Paneth-like": "Paneth cells",
                "Stem Cell": "Stem Cell",
                "TA": "TA",
                "Enteriendocrine": "Enteroendocrine cells",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "ileum", "wang20_ileum.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["n_counts"].values[:, None]))\
                                   .multiply(1 / 10000)
