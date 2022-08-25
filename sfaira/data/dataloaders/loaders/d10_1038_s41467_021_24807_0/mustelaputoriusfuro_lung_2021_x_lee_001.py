import os

import pandas as pd

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, "GSE171828_Seurat_object_total_cells.Rds.gz")
    # Fine grained annotation of myleoid:
    fn_meta = os.path.join(data_dir, "GSE171828_monocyte_macrophage_DC_annotation.tsv.gz")
    fn = buffered_decompress(fn)
    fn_meta = buffered_decompress(fn_meta)
    anndata2ri.activate()
    adata = r(
        f"library(Seurat)\n"
        f"obj <- readRDS('{fn}')\n"
        f"new_obj = CreateSeuratObject(counts = obj@assays$RNA@counts)\n"
        f"names(rownames(new_obj@assays$RNA@counts)) <- NULL\n"
        f"new_obj@meta.data = obj@meta.data\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    tab_meta = pd.read_csv(fn_meta, sep="\t", index_col=0)
    adata.obs["Annotation_fine"] = [
        tab_meta.loc[x, "Annotation"] if x in tab_meta.index else y
        for x, y in zip(adata.obs_names, adata.obs["Annotation"].values)]
    return adata
