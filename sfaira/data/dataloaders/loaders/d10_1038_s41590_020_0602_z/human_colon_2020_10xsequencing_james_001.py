import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/james20.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc = "10x technology"
        self.author = "James"
        self.disease = "healthy"
        self.doi = "10.1038/s41590-020-0602-z"
        self.normalization = "raw"
        self.organ = "colon"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2020

        self.gene_id_symbols_var_key = "index"
        self.gene_id_ensembl_var_key = "gene_ids"
        self.cell_types_original_obs_key = "cell_type"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "james20.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

    return adata
