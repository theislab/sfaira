import anndata
import os
import pandas as pd


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    fn_meta = os.path.join(data_dir, "complete_MetaData_70170cells_scp.csv")
    adata = anndata.read_csv(fn).T
    tab_meta = pd.read_csv(fn_meta, index_col=0, header=0, skiprows=[1])
    # Some cell identifiers occur repeatedly in meta data file, because they are listed for different figures:
    del tab_meta["dataset"]
    tab_meta = tab_meta.drop_duplicates()
    # A few cells are not be annotated, we discard these here.
    adata = adata[[x in tab_meta.index for x in adata.obs_names], :].copy()
    adata.obs = tab_meta.loc[adata.obs_names, :]

    return adata
