import anndata
import os

import pandas as pd


def load(data_dir, sample_fn, **kwargs):
    fn_x = os.path.join(data_dir, "SCP975", "expression", sample_fn)
    fn_meta = os.path.join(data_dir, "SCP975", "cluster", sample_fn.replace("pre_EF", "tSNE"))
    tab_meta = pd.read_csv(fn_meta, index_col=0)
    tab_meta = tab_meta.iloc[1:, :].copy()
    adata = anndata.read_csv(fn_x).T
    adata = adata[tab_meta.index, :].copy()
    adata.obs = tab_meta
    return adata
