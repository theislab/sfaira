import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "lukassen20_lung_orig.processed.h5ad",
    "lukassen20_airway_orig.processed.h5ad"
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
        self.id = f"human_lung_2020_10x_lukassen_{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_" \
                  f"10.1101/2020.03.13.991455"

        self.download_url_data = f"https://covid19.cog.sanger.ac.uk/{self.sample_fn}"
        self.download_url_meta = None

        self.author = "Lukassen"
        self.doi = "10.1101/2020.03.13.991455"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

    def _load(self):
        fn = os.path.join(self.data_dir, self.sample_fn)
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["nCount_RNA"].values[:, None])).multiply(1 / 10000)
        self.set_unknown_class_id(ids=["1_Unicorns and artifacts"])

        return adata
