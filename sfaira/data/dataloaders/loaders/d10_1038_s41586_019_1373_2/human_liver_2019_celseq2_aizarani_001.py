import anndata
import os
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5FNormalhumanlivercellatlasdata%2Etxt%2Egz"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5Fclusterpartition%2Etxt%2Egz"

        self.assay_sc = "CEL-seq2"
        self.author = "Aizarani"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41586-019-1373-2"
        self.normalization = "raw"
        self.sample_source = "primary_tissue"
        self.organ = "liver"
        self.organism = "human"
        self.year = 2019

        self.gene_id_symbols_var_key = "index"
        self.cell_type_obs_key = "CellType"

        self.set_dataset_id(idx=1)


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
