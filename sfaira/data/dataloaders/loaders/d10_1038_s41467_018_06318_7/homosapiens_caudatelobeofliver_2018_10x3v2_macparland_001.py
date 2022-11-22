import anndata
import os
import pandas as pd


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE115469.csv.gz"),
        os.path.join(data_dir, "GSE115469_labels.txt")
    ]
    adata = anndata.read_csv(fn[0]).T
    celltype_df = pd.read_csv(fn[1], sep="\t").set_index("CellName")
    adata.obs["celltype"] = [str(celltype_df.loc[i]["Cluster#"]) for i in adata.obs.index]

    return adata
