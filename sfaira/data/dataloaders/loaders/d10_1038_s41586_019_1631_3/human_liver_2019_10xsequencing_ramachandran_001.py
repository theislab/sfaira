import os
import numpy as np


def load(data_dir, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, "tissue.rdata")
    anndata2ri.activate()  # TODO: remove global activation of anndata2ri and use localconverter once it's fixed
    adata = r(
        f"library(Seurat)\n"
        f"load('{fn}')\n"
        f"new_obj = CreateSeuratObject(counts = tissue@raw.data)\n"
        f"new_obj@meta.data = tissue@meta.data\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    adata.obs["nGene"] = adata.obs["nGene"].astype(np.int32)

    return adata
