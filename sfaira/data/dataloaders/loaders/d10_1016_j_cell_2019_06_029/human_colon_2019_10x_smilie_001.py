import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_colon_2019_10x_smilie_001_10.1016/j.cell.2019.06.029"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/smillie19_epi.processed.h5ad"
        self.download_url_meta = None

        self.author = "Smilie"
        self.doi = "10.1016/j.cell.2019.06.029"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "colonic epithelium"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

    def _load(self):
        fn = os.path.join(self.data_dir, "smillie19_epi.processed.h5ad")
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

        return adata
