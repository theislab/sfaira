import anndata
import os
import scipy.sparse


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = scipy.sparse.csc_matrix(adata.X)

    return adata
