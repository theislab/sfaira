import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/lukowski19.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Lukowski"
        self.disease = "healthy"
        self.doi_journal = "10.15252/embj.2018100811"
        self.normalization = "raw"
        self.organ = "retina"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.year = 2019

        self.gene_id_symbols_var_key = "index"
        self.gene_id_ensembl_var_key = "gene_ids"
        self.cell_type_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "lukowski19.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

    return adata
