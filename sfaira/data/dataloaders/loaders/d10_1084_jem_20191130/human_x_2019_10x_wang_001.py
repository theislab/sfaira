import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "wang20_colon.processed.h5ad",
    "wang20_ileum.processed.h5ad",
    "wang20_rectum.processed.h5ad"
]


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def __init__(
            self,
            sample_fn: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, sample_fns=SAMPLE_FNS, data_path=data_path, meta_path=meta_path,
                         cache_path=cache_path, **kwargs)
        self.download_url_data = f"https://covid19.cog.sanger.ac.uk/{self.sample_fn}"
        self.download_url_meta = None

        organ = self.sample_fn.split("_")[1].split(".")[0]

        self.author = "Wang"
        self.doi = "10.1084/jem.20191130"
        self.healthy = True
        self.normalization = "raw"
        self.organ = organ
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        if organ == "colon":
            self.class_maps = {
                "0": {
                    "Progenitor": "Enterocyte Progenitors",
                    "Enterocyte": "Enterocytes",
                    "Goblet": "Goblet cells",
                    "TA": "TA",
                    "Paneth-like": "Paneth cells",
                    "Stem Cell": "Stem cells",
                    "Enteriendocrine": "Enteroendocrine cells",
                },
            }
        elif organ == "ileum":
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
        else:
            self.class_maps = {
                "0": {
                    "Progenitor": "Enterocyte progenitor",
                    "Goblet": "Goblet",
                    "Enterocyte": "Enterocyte",
                    "Paneth-like": "Paneth-like",
                    "Stem Cell": "Stem Cell",
                    "TA": "TA",
                    "Enteriendocrine": "Enteroendocrine",
                },
            }

    def _load(self):
        fn = os.path.join(self.data_dir, self.sample_fn)
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

        return adata
