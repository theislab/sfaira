import anndata
import os


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "Full_obj_raw_counts_nosoupx_v2.h5ad")
    adata = anndata.read_h5ad(fn)
    return adata
