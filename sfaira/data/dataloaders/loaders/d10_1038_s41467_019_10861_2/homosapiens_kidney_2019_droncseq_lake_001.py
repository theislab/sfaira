import anndata
import os
import pandas as pd


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv.gz"),
        os.path.join(data_dir, "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotations.csv.gz")
    ]
    adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t").T)
    annot = pd.read_csv(fn[1], index_col=0, dtype="category")
    adata.obs["celltype"] = [annot.loc[i.split("_")[0][1:]]["Annotation"] for i in adata.obs.index]

    return adata
