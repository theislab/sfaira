import anndata
import os
import scipy.sparse


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "my-data-file.h5ad")
    adata = anndata.read(fn)
    adata.X = scipy.sparse.csr_matrix(adata.X)

    return adata
