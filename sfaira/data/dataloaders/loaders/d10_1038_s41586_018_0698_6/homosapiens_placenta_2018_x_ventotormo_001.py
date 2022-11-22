import os
import pandas as pd
import anndata


def load(data_dir, sample_fn, **kwargs):
    fn = [
        os.path.join(data_dir, f"{sample_fn}.1.zip"),
        os.path.join(data_dir, f"{sample_fn}.2.zip"),
    ]
    adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t", index_col="Gene").T)
    df = pd.read_csv(fn[1], sep="\t")
    for i in df.columns:
        adata.obs[i] = [df.loc[j][i] for j in adata.obs.index]

    adata.var["ensembl"] = [i.split("_")[1] for i in adata.var.index]
    adata.var["names"] = [i.split("_")[0] for i in adata.var.index]
    adata.var = adata.var.reset_index().reset_index().drop("index", axis=1)
    adata = adata[:, ~adata.var.index.isin(
        ["", "-1", "-10", "-11", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "A.2", "A.3"])].copy()

    return adata
