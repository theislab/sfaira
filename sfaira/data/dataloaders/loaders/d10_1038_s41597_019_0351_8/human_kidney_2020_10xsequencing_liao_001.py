import anndata
import os
import pandas as pd
import scipy.io
import gzip
import tarfile

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131685/suppl/GSE131685_RAW.tar"
        self.download_url_meta = None

        self.assay_sc = "10x technology"
        self.author = "Liao"
        self.disease = "healthy"
        self.normalization = "raw"
        self.organ = "kidney"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.year = 2020
        self.doi = "10.1038/s41597-019-0351-8"

        self.gene_id_symbols_var_key = "names"
        self.gene_id_ensembl_var_key = "ensembl"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE131685_RAW.tar")
    adatas = []
    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            if "_matrix.mtx.gz" in member.name:
                name = "_".join(member.name.split("_")[:-1])
                with gzip.open(tar.extractfile(member), "rb") as mm:
                    X = scipy.io.mmread(mm).T.tocsr()
                obs = pd.read_csv(tar.extractfile(name + "_barcodes.tsv.gz"), compression="gzip", header=None,
                                  sep="\t", index_col=0)
                obs.index.name = None
                var = pd.read_csv(tar.extractfile(name + "_features.tsv.gz"), compression="gzip", header=None,
                                  sep="\t").iloc[:, :2]
                var.columns = ["ensembl", "names"]
                var.index = var["ensembl"].values
                adata = anndata.AnnData(X=X, obs=obs, var=var)
                adata.obs["sample"] = name
                adatas.append(adata)
    adata = adatas[0].concatenate(adatas[1:])

    return adata
