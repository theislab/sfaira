import anndata
import os
import pandas as pd


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE124395_Normalhumanlivercellatlasdata.txt.gz"),
        os.path.join(data_dir, "GSE124395_clusterpartition.txt.gz")
    ]
    adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t").T)
    celltype_df = pd.read_csv(fn[1], sep=" ")
    adata = adata[[i in celltype_df.index for i in adata.obs.index]].copy()
    adata.obs["CellType"] = [str(celltype_df.loc[i]["sct@cpart"]) for i in adata.obs.index]

    return adata
