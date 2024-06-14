import anndata
import os
import scipy.sparse

def load(data_dir,sample_fn, **kwargs):
    import anndata2ri
    from rpy2.robjects import r
    anndata2ri.activate()
    fn = os.path.join(data_dir, sample_fn)
    adata = r(
        f"library(Seurat)\n"
        f"new_obj = readRDS('{fn}')))\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    return adata
