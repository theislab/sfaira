import anndata
import os
import pandas as pd


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE114397_HMLE_TGFb.tsv.gz")
    tab = pd.read_csv(fn, sep="\t")
    adata = anndata.AnnData(
        tab.iloc[:, 1:].values,
        obs=pd.DataFrame({"sample": tab.iloc[:, 0].values}),
        var=pd.DataFrame({}, index=tab.columns[1:]))
    return adata
