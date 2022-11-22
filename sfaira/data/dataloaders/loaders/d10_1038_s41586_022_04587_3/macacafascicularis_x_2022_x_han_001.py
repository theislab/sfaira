import anndata
import os
import pandas as pd
import scipy.sparse


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    fn_meta = os.path.join(data_dir, "Metadata_ALL_Tissue_Global_Clustering.txt.gz")
    tab_x = pd.read_csv(fn, header=0, index_col=0, compression="gzip", sep="\t").T
    tab_meta = pd.read_csv(fn_meta, header=0, index_col=False, compression="gzip", sep="\t")
    tab_meta.index = tab_meta["ID"].values
    adata = anndata.AnnData(scipy.sparse.csr_matrix(tab_x.values), var=pd.DataFrame({}, index=tab_x.columns))
    adata.obs_names = tab_x.index
    adata.obs = tab_meta.loc[adata.obs_names, :]
    return adata
