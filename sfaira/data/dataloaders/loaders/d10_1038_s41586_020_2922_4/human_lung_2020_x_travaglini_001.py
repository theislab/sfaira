import anndata
import os
import scipy.sparse
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    if sample_fn.split("_")[0] == "droplet":
        norm_const = 10000
        sf_key = "nUMI"
    else:
        norm_const = 1000000
        sf_key = "nReads"
    adata = anndata.read(fn)
    adata.X = scipy.sparse.csc_matrix(adata.X)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs[sf_key].values[:, None])).multiply(1 / norm_const)

    return adata
