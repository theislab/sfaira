import anndata
import os
import pandas as pd


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE135893_matrix.mtx.gz"),
        os.path.join(data_dir, "GSE135893_genes.tsv.gz"),
        os.path.join(data_dir, "GSE135893_barcodes.tsv.gz"),
        os.path.join(data_dir, "GSE135893_IPF_metadata.csv.gz"),
        os.path.join(data_dir, "aba1972_table_s2.csv"),
    ]
    adata = anndata.read_mtx(fn[0]).T
    adata.var = pd.read_csv(fn[1], index_col=0, header=None, names=["ids"])
    adata.obs = pd.read_csv(fn[2], index_col=0, header=None, names=["barcodes"])
    obs = pd.read_csv(fn[3], index_col=0)
    obs2 = pd.read_csv(fn[4], index_col=0)
    obs["Chemistry"] = [obs2.loc[x, "Chemistry"] for x in obs["orig.ident"].values]
    obs["Gender"] = [obs2.loc[x, "Gender"] for x in obs["orig.ident"].values]
    adata = adata[obs.index.tolist(), :].copy()
    adata.obs = obs

    return adata
