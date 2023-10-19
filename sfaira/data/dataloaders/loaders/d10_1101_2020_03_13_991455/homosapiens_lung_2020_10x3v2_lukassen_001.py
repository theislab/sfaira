import anndata
import os
import numpy as np
import scipy.sparse


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["nCount_RNA"].values[:, None])).multiply(1 / 10000)

    return adata
