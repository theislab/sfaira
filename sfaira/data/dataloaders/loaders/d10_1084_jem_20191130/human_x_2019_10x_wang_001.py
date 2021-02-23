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
        super().__init__(sample_fn=sample_fn, data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        organ = self.sample_fn.split("_")[1].split(".")[0]
        self.id = f"human_{organ}_2019_10x_wang_{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_10.1084/jem.20191130"

        self.download_url_data = f"https://covid19.cog.sanger.ac.uk/wang20_{organ}.processed.h5ad"
        self.download_url_meta = None

        self.author = "Wang"
        self.doi = "10.1084/jem.20191130"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "colon" if organ == "colon" else "ileum" if organ == "ileum" else "rectum"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

    def _load(self):
        fn = os.path.join(self.data_dir, self.sample_fn)
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

        return adata
