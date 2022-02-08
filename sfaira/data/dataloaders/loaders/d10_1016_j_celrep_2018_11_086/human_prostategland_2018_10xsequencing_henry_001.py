import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    ToDo: revisit these cell type maps, Club and Hillock are described in this paper.
      Club,epithelial cell of prostate
      Hillock,epithelial cell of prostate
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/henry18_0.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Henry"
        self.disease = "healthy"
        self.doi_journal = "10.1016/j.celrep.2018.11.086"
        self.layer_counts = "X"
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.organ = "prostate gland"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.year = 2018

        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"

        self.cell_type_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "henry18_0.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

    return adata
