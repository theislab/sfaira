import anndata
import os
import numpy as np
import scipy.sparse


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "miller20.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["nUMI"].values[:, None])).multiply(1 / 10000)

    return adata
