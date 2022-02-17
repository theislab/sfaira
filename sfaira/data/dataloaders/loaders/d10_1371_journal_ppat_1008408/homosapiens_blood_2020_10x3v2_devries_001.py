import os
import pandas as pd
import numpy as np

def load(data_dir, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, "seurat_object_final_anonymized.Rds")
    anndata2ri.activate()
    adata = r(
        f"library(Seurat)\n"
        f"so = readRDS('{fn}')\n"
        f"so_updated = UpdateSeuratObject(so)\n"
        f"as.SingleCellExperiment(so_updated)\n"
    )
    adata.obs["stimulation"] = pd.Categorical([i if i == "unstimulated" else "C. albicans" for i in adata.obs["stimulation"]])
    adata.obs["nFeature_RNA"] = adata.obs["nFeature_RNA"].astype(np.int32)
    return adata
