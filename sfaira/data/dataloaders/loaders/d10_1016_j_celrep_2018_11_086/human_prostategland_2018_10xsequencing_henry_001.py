import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    ToDo: revisit these cell type maps, Club and Hillock are described in this paper.
      Club,epithelial cell of prostate
      Hillock,epithelial cell of prostate
    """

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/henry18_0.processed.h5ad"
        self.download_url_meta = None

        self.author = "Henry"
        self.doi = "10.1016/j.celrep.2018.11.086"
        self.healthy = True
        self.normalization = "raw"
        self.state_exact = "healthy"
        self.organ = "prostate gland"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.year = 2018

        self.var_symbol_col = "index"
        self.obs_key_cellontology_original = "CellType"

        self.set_dataset_id(idx=1)

    def _load(self):
        fn = os.path.join(self.data_dir, "henry18_0.processed.h5ad")
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

        return adata
