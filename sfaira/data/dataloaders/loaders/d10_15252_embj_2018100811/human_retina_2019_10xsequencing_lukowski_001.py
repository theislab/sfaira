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

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/lukowski19.processed.h5ad"
        self.download_url_meta = None

        self.author = "Lukowski"
        self.doi = "10.15252/embj.2018100811"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "retina"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"
        self.var_ensembl_col = "gene_ids"
        self.obs_key_cellontology_original = "CellType"

        self.set_dataset_id(idx=1)

        self.set_unknown_class_id(ids=["unannotated", "unspecified"])


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "lukowski19.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

    return adata
