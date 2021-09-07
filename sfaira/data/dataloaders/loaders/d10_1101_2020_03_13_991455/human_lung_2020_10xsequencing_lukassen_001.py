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

        self.assay_sc = "10x 3' v2"
        self.author = "Lukassen"
        self.disease = "healthy"
        self.doi_journal = "10.15252/embj.20105114"
        self.doi_preprint = "10.1101/2020.03.13.991455"
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.year = 2020

        self.gene_id_symbols_var_key = "index"
        self.cell_type_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["nCount_RNA"].values[:, None])).multiply(1 / 10000)

    return adata
