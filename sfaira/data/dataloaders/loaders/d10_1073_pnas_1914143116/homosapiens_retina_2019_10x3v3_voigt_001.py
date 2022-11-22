import anndata
import os
import numpy as np


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "voigt19.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)

    return adata
