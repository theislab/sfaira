import anndata
import os
import scipy.sparse


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.obs["age"] = [f"{x}-year-old human stage" for x in adata.obs["age"].values]

    return adata
