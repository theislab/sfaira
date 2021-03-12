import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad"
        self.download_url_meta = None

        self.author = "Habib"
        self.doi = "10.1038/nmeth.4407"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "brain"
        self.organism = "human"
        self.assay_sc = "DroNc-seq"
        self.state_exact = "healthy"
        self.year = 2017
        self.sample_source = "primary_tissue"

        self.var_symbol_col = "index"
        self.cellontology_original_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "habib17.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

    return adata
