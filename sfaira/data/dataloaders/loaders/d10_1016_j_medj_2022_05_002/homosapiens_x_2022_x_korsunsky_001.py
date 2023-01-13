import anndata
import os
import scipy.io

import pandas as pd


def load(data_dir, **kwargs):
    fn_obs = os.path.join(data_dir, "umis_barcodes.tsv")
    fn_var = os.path.join(data_dir, "umis_genes.tsv")
    fn_x = os.path.join(data_dir, "umis.mtx")
    fn_meta = os.path.join(data_dir, "metaData.txt")
    x = scipy.io.mmread(fn_x).T.tocsr()
    tab_obs = pd.read_csv(fn_obs, index_col=0, header=None, sep="\t", dtype=str)
    tab_var = pd.read_csv(fn_var, index_col=0, header=None, sep="\t", dtype=str)
    tab_var.index = [str(x) for x in tab_var.index]
    tab_meta = pd.read_csv(fn_meta, index_col=0, header=0, skiprows=[1], sep="\t", dtype=str)
    tab_meta["cell_type_within_tissue"] = [str(x) for x in tab_meta["cell_type_within_tissue"].values]
    adata = anndata.AnnData(x, obs=tab_obs, var=tab_var)
    adata.obs = tab_meta
    return adata
