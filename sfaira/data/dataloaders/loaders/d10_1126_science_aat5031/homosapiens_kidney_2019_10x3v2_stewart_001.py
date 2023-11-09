import anndata
import os
import numpy as np


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "Mature_Full_v2.1.h5ad"),
        os.path.join(data_dir, "Fetal_full.h5ad")
    ]
    adult = anndata.read_h5ad(fn[0])
    fetal = anndata.read_h5ad(fn[1])
    # TODO this is is not a controlled field
    adult.obs["development"] = "adult"
    fetal.obs["development"] = "fetal"
    adata = adult.concatenate(fetal)
    adata.X = np.expm1(adata.X)

    return adata
