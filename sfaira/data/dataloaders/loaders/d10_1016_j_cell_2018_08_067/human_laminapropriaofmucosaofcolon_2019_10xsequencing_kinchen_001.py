import os
import pandas as pd
import anndata as ad
import scipy.sparse
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    fn = [
        os.path.join(data_dir, f"GSE114374_Human_{sample_fn}_expression_matrix.txt.gz"),
        os.path.join(data_dir, f"{sample_fn.lower()}_meta_data_stromal_with_donor.txt"),
    ]
    matrix = pd.read_csv(fn[0], sep="\t")
    obs = pd.read_csv(fn[1], sep="\t", index_col=3)
    adata = ad.AnnData(matrix.T)
    adata.X = scipy.sparse.csc_matrix(np.expm1(adata.X))
    adata.obs = obs
    s_dict = {"F": "female", "M": "male"}
    adata.obs['Sex'] = [s_dict[i] for i in adata.obs['Sex']]
    adata.obs['Age'] = [str(x) + "-year-old homosapiens stage" for x in adata.obs['Age'].values]

    return adata
