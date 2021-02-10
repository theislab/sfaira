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
        self.id = "human_testis_2018_10x_guo_001_10.1038/s41422-018-0099-2"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/guo18_donor.processed.h5ad"
        self.download_url_meta = None

        self.author = "Cairns"
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

        self.class_maps = {
            "0": {
                "Elongated Spermatids": "Elongated Spermatids",
                "Leydig cells": "Leydig cells",
                "Early Primary Spermatocytes": "Early Primary Spermatocytes",
                "Round Spermatids": "Round Spermatids",
                "Endothelial cells": "Endothelial cells",
                "Macrophages": "Macrophages",
                "Myoid cells": "Myoid cells",
                "Differentiating Spermatogonia": "Differentiating Spermatogonia",
                "Late primary Spermatocytes": "Late primary Spermatocytes",
                "Spermatogonial Stem cell": "Spermatogonial Stem cell",
                "Sertoli cells": "Sertoli cells",
            },
        }

    def _load(self):
        fn = os.path.join(self.doi_path, "guo18_donor.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["n_counts"].values[:, None]))\
                                   .multiply(1 / 10000)
