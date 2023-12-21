import anndata
import os
import scipy.sparse
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read_h5ad(fn)
    adata.X = np.round((np.expm1(adata.X.A) * adata.obs["nCount_RNA"][:, None] / 10000))
    adata.X = scipy.sparse.csr_matrix(adata.X)

    return adata
