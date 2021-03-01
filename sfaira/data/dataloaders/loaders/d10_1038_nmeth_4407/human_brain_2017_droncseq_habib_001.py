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
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad"
        self.download_url_meta = None

        self.author = "Habib"
        self.doi = "10.1038/nmeth.4407"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "brain"
        self.organism = "human"
        self.protocol = "DroNc-seq"
        self.state_exact = "healthy"
        self.year = 2017

        self.var_symbol_col = "index"
        self.obs_key_cellontology_original = "CellType"

        self.set_dataset_id(idx=1)

        self.set_unknown_class_id(ids=["Unclassified"])


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "habib17.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

    return adata
