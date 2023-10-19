import anndata
import os
import pandas as pd
import scipy.io

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, **kwargs):
    data_dir = buffered_decompress(os.path.join(data_dir, "1M_PBMC_T1D_Parse.zip"))
    fn_x = os.path.join(data_dir, "DGE_1M_PBMC.mtx")
    fn_obs = os.path.join(data_dir, "cell_metadata_1M_PBMC.csv")
    fn_var = os.path.join(data_dir, "all_genes_1M_PBMC.csv")
    obs = pd.read_csv(fn_obs, index_col=False)
    var = pd.read_csv(fn_var, index_col=0)
    x = scipy.io.mmread(fn_x).tocsr()
    adata = anndata.AnnData(x, obs=obs, var=var)
    return adata
