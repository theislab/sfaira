import anndata
import os
import pandas as pd


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "E-MTAB-5061.processed.1.zip"),
        os.path.join(data_dir, "E-MTAB-5061.sdrf.txt")
    ]
    df = pd.read_csv(fn[0], sep="\t")
    df.index = df.index.get_level_values(0)
    df = df.drop("#samples", axis=1)
    df = df.T.iloc[:, :26178]
    adata = anndata.AnnData(df)
    adata.obs = pd.read_csv(fn[1], sep="\t").set_index("Source Name").loc[adata.obs.index]
    # filter observations which are not cells (empty wells, low quality cells etc.)
    adata = adata[adata.obs["Characteristics[cell type]"] != "not applicable"].copy()

    return adata
