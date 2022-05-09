import anndata
import numpy as np
import pandas as pd
import os

from scipy.sparse import csr_matrix


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    # reconstruct raw counts
    adata.layers['counts'] = adata.X.expm1().multiply(csr_matrix(adata.obs['size_factors'].to_numpy().reshape((-1, 1))))
    adata.obs["age"] = pd.Categorical(["nan" if np.isnan(x) else str(x) for x in adata.obs["age"].values])

    return adata
