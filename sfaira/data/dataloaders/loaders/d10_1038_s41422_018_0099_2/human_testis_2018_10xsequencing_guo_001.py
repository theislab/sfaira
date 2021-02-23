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
        self.id = "human_testis_2018_10xsequencing_guo_001_10.1038/s41422-018-0099-2"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/guo18_donor.processed.h5ad"
        self.download_url_meta = None

        self.author = "Guo"
        self.doi = "10.1038/s41422-018-0099-2"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "testis"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

    def _load(self):
        fn = os.path.join(self.data_dir, "guo18_donor.processed.h5ad")
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

        return adata
