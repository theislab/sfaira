import anndata
import os
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    TODO:  add disease from status and diagnosis fields, healthy is "control"
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fmatrix%2Emtx%2Egz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fgenes%2Etsv%2Egz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fbarcodes%2Etsv%2Egz"
        ]
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5FIPF%5Fmetadata%2Ecsv%2Egz"

        self.author = "Habermann"
        self.doi = "10.1101/753806"
        self.normalization = "raw"
        self.organ = "lung parenchyma"
        self.organism = "human"
        self.assay_sc = "10x technology"
        self.year = 2020
        self.sample_source = "primary_tissue"

        self.gene_id_symbols_var_key = "index"

        self.cell_types_original_obs_key = "celltype"
        self.state_exact_obs_key = "Diagnosis"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE135893_matrix.mtx.gz"),
        os.path.join(data_dir, "GSE135893_genes.tsv.gz"),
        os.path.join(data_dir, "GSE135893_barcodes.tsv.gz"),
        os.path.join(data_dir, "GSE135893_IPF_metadata.csv.gz"),
    ]
    adata = anndata.read_mtx(fn[0]).T
    adata.var = pd.read_csv(fn[1], index_col=0, header=None, names=["ids"])
    adata.obs = pd.read_csv(fn[2], index_col=0, header=None, names=["barcodes"])
    obs = pd.read_csv(fn[3], index_col=0)
    adata = adata[obs.index.tolist(), :].copy()
    adata.obs = obs

    return adata
