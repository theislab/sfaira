import os


def load(data_dir, sample_fn, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, "covid19.macrophages", sample_fn)
    object_name = sample_fn.replace(".Rds", "")
    anndata2ri.activate()
    adata = r(
        f"library(Seurat)\n"
        f"readRDS('{fn}')\n"
        f"new_obj = CreateSeuratObject(counts = `{object_name}`@assays$RNA@counts)\n"
        f"new_obj@meta.data = `{object_name}`@meta.data\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    # TODO SarsCov2 counts are not yet loaded here, they are in
    # `{object_name}`@assays$SCoV2

    return adata
