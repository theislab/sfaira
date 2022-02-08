import anndata
import os
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121862/suppl/" \
            "GSE121862%5FUCSD%2DWU%5FSingle%5FNuclei%5FCluster%5FAnnotated%5FRaw%5FUMI%5FMatrix%2Etsv%2Egz"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121862/suppl/" \
            "GSE121862%5FUCSD%2DWU%5FSingle%5FNuclei%5FCluster%5FAnnotations%2Ecsv%2Egz"

        self.assay_sc = "DroNc-seq"
        self.author = "Lake"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41467-019-10861-2"
        self.layer_counts = "X"
        self.organ = "kidney"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.year = 2019

        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"
        self.cell_type_obs_key = "celltype"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv.gz"),
        os.path.join(data_dir, "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotations.csv.gz")
    ]
    adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t").T)
    annot = pd.read_csv(fn[1], index_col=0, dtype="category")
    adata.obs["celltype"] = [annot.loc[i.split("_")[0][1:]]["Annotation"] for i in adata.obs.index]

    return adata
