import anndata
import os


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "menon19.processed.h5ad")
    adata = anndata.read(fn)

    return adata
