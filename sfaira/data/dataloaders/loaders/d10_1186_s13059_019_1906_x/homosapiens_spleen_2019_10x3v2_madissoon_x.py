import anndata
import os
import scipy.sparse


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    if sample_fn != "madissoon19_lung.processed.h5ad":
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None]))\
            .multiply(1 / 10000)
    # Cell type column called differently in madissoon19_lung.processed.h5ad:
    if sample_fn == "madissoon19_lung.processed.h5ad":
        adata.obs["Celltypes"] = adata.obs["CellType"]
        del adata.obs["CellType"]

    return adata
