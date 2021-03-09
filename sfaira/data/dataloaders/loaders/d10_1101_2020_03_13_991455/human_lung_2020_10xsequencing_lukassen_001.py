import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase

SAMPLE_FNS = [
    "lukassen20_lung_orig.processed.h5ad",
    "lukassen20_airway_orig.processed.h5ad"
]


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = f"https://covid19.cog.sanger.ac.uk/{self.sample_fn}"
        self.download_url_meta = None

        self.author = "Lukassen"
        self.doi = "10.1101/2020.03.13.991455"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.assay_sc = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2020
        self.sample_source = "primary_tissue"

        self.var_symbol_col = "index"
        self.cellontology_original_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["nCount_RNA"].values[:, None])).multiply(1 / 10000)

    return adata
