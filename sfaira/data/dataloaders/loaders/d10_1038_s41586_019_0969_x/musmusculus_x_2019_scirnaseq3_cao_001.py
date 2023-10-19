import anndata
import os
import pandas as pd
import scipy.io

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, **kwargs):
    fn_x = os.path.join(data_dir, "GSE119945_gene_count.txt.gz")  # Note: This is an .mtx!
    fn_x = buffered_decompress(fn_x)
    fn_obs = os.path.join(data_dir, "GSE119945_cell_annotate.csv.gz")
    fn_var = os.path.join(data_dir, "GSE119945_gene_annotate.csv.gz")
    fn_obs2 = os.path.join(data_dir, "cell_annotate.csv")
    # replace the simple data loading code below with the code required to load your data file(s)
    tab_obs = pd.read_csv(fn_obs, index_col=0, sep=",", compression="gzip", dtype=str)
    tab_obs2 = pd.read_csv(fn_obs2, index_col=0, sep=",", dtype=str)
    tab_obs2["Main_cell_type"] = [str(x) for x in tab_obs2["Main_cell_type"].values]
    tab_obs2["detected_doublet"] = [str(x) for x in tab_obs2["detected_doublet"].values]
    tab_obs2["development_stage"] = [str(x) for x in tab_obs2["development_stage"].values]
    tab_obs2["embryo_sex"] = [str(x) for x in tab_obs2["embryo_sex"].values]
    tab_var = pd.read_csv(fn_var, index_col=0, sep=",", compression="gzip")
    x = scipy.io.mmread(fn_x).T.tocsr()
    adata = anndata.AnnData(x, var=tab_var)
    adata.obs_names = tab_obs.index
    adata.obs = tab_obs2.loc[adata.obs_names, :]
    return adata
