import anndata
import os


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read_h5ad(fn)

    return adata
