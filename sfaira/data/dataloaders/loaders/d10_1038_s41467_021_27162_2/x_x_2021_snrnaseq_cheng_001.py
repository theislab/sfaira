import os


def load(data_dir, sample_fn, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, sample_fn)
    anndata2ri.activate()
    adata = r(
        f"library(Seurat)\n"
        f"obj <- readRDS('{fn}')\n"
        f"new_obj = CreateSeuratObject(counts = obj@assays$RNA@counts)\n"
        f"new_obj@meta.data = obj@meta.data\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    print(adata.obs.head())

    return adata
