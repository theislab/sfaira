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
        self.id = "human_lung_2020_10x_miller_001_10.1016/j.devcel.2020.01.033"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/miller20.processed.h5ad"
        self.download_url_meta = None

        self.author = "Miller"
        self.doi = "10.1016/j.devcel.2020.01.033"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "Cell_type"

    def _load(self):
        fn = os.path.join(self.data_dir, "miller20.processed.h5ad")
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["nUMI"].values[:, None])).multiply(1 / 10000)
        self.set_unknown_class_id(ids=[
            "Bud tip adjacent",
            "Bud tip progenitor",
            "Submucosal gland",
            "Submucosal gland basal",
        ])

        return adata
