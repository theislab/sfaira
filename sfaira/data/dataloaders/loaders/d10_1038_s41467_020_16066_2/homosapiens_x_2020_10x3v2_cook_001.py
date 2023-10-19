import anndata
import os

import pandas as pd


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    sample_fn_prefix = sample_fn.replace("UMI_matrix.csv.gz", "")
    meta_fn_suffix = "metadata.csv.gz"
    fn_meta = os.path.join(data_dir, sample_fn_prefix + meta_fn_suffix)
    adata = anndata.read_csv(fn).T
    tab = pd.read_csv(fn_meta, index_col=0)
    # In some datasets, extra coordinates are in the data file, these are discarded here:
    adata = adata[[x in tab.index for x in adata.obs_names], :].copy()
    adata.obs = tab.loc[adata.obs_names, :]

    return adata
