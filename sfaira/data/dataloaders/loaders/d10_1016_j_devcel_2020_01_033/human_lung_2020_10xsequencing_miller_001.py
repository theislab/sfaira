import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/miller20.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Miller"
        self.disease = "healthy"
        self.doi_journal = "10.1016/j.devcel.2020.01.033"
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2020

        self.gene_id_symbols_var_key = "index"
        self.cell_type_obs_key = "Cell_type"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "miller20.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["nUMI"].values[:, None])).multiply(1 / 10000)

    return adata
