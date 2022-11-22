import anndata
import os


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom")
    adata = anndata.read_loom(fn)
    adata = adata[adata.obs["emptydrops_is_cell"] == "t"].copy()

    return adata
