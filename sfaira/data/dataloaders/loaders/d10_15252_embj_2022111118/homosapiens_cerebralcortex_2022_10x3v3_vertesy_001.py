import os
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, "GSE205554_integration_all_preGruffi.Rds.gz")
    anndata2ri.activate()
    adata = r(
        f"library(Seurat)\n"
        f"as.SingleCellExperiment(readRDS('{fn}'), assay='RNA')\n"
    )
    adata.X = adata.X.astype(np.float32)
    adata.layers["logcounts"] = adata.layers["logcounts"].astype(np.float32)
    adata = adata[adata.obs["project"] == "Ver"].copy()
    adata.obs["organoid_age_days"] = "120"

    return adata
