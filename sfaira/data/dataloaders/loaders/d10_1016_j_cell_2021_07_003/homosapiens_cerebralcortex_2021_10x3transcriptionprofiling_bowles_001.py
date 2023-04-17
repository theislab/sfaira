import os
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, f"All_{sample_fn}_Clustered_Annotated.rds")
    anndata2ri.activate()
    adata = r(
        f"library(Seurat)\n"
        f"so = readRDS('{fn}')\n"
        f"as.SingleCellExperiment(so, assay='RNA')\n"
    )
    adata.X = adata.X.astype(np.float32)
    age_dict = {
        "2mo": "60",
        "4mo": "120",
        "6mo": "180",
    }
    adata.obs["organoid_age_days"] = age_dict[sample_fn]
    adata.obs["Run"] = adata.obs["Run"].astype(str)

    return adata
