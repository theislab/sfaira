import anndata
import os
# import any packages you require for dataloading here. you can assume packages like numpy and pandas being available


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "fetal_liver_alladata_.h5ad")
    adata = anndata.read_h5ad(fn)

    return adata
