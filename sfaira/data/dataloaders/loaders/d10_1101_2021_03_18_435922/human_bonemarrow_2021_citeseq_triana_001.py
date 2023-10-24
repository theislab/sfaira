import os


def load(data_dir, sample_fn, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    object_name = "Integration_200AB" if sample_fn == "200AB" else sample_fn
    fn = os.path.join(data_dir, sample_fn + ".rda")
    anndata2ri.activate()  # TODO: remove global activation of anndata2ri and use localconverter once it's fixed
    # file_name is already a Seurat object, maybe not optimally solved here, need a way to select correct count matrix.
    adata = r(
        f"library(Seurat)\n"
        f"load('{fn}')\n"
        f"new_obj = CreateSeuratObject(counts = {object_name}@assays$RNA@counts)\n"
        f"new_obj@meta.data = {object_name}@meta.data\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    if "class" not in adata.obs_keys():
        adata.obs["class"] = "healthy"
    return adata
