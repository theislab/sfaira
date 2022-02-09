import anndata
import os
import pandas as pd
import scipy.io

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = [
            'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_scRNAseq_raw_UMIs.mtx',
            'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_scRNAseq_genes.tsv',
            'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_scRNAseq_barcodes.tsv']
        self.download_url_meta = 'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_metadata.txt'

        self.assay_sc = "10x 3' v2"
        self.author = "Cano-Gamez"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41467-020-15543-y"
        self.layer_counts = "X"
        self.organ = "blood"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.year = 2020

        self.feature_symbol_var_key = "names"
        self.feature_type = "rna"
        self.cell_type_obs_key = "cluster.id"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [os.path.join(data_dir, "NCOMMS-19-7936188_scRNAseq_raw_UMIs.mtx"),
          os.path.join(data_dir, 'NCOMMS-19-7936188_scRNAseq_barcodes.tsv'),
          os.path.join(data_dir, 'NCOMMS-19-7936188_scRNAseq_genes.tsv'),
          os.path.join(data_dir, "NCOMMS-19-7936188_metadata.txt")]

    with open(fn[0], 'rb') as mm:
        X = scipy.io.mmread(mm).T.tocsr()
    obs = pd.read_csv(fn[1], header=None, sep="\t", index_col=0)
    obs.index.name = None
    var = pd.read_csv(fn[2], header=None, sep="\t", names=['names'])
    var.index = var['names'].values

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    metadata = pd.read_csv(fn[3], sep="\t")
    adata.obs = adata.obs.join(metadata, how="inner")
    return adata
