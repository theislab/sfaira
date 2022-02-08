import anndata
import os
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fraw%5Fcounts%2Ecsv%2Egz"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fbarcodes%5Fcell%5Ftypes%2Etxt%2Egz"

        self.assay_sc = "Drop-seq"
        self.author = "Braga"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41591-019-0468-5"
        self.layer_counts = "X"
        self.organ = "lung"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.state_exact = "uninvolved areas of tumour resection material"
        self.year = 2019

        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"
        self.cell_type_obs_key = "celltype"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE130148_raw_counts.csv.gz"),
        os.path.join(data_dir, "GSE130148_barcodes_cell_types.txt.gz"),
    ]
    adata = anndata.read_csv(fn[0]).T
    adata.obs = pd.read_csv(fn[1], sep="\t", index_col=0)

    return adata
