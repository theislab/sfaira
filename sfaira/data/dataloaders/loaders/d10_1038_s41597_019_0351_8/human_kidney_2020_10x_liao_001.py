import anndata
import os
from typing import Union
import pandas as pd
import scipy.io
import gzip
import tarfile

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_kidney_2020_10x_liao_001_10.1038/s41597-019-0351-8"

        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131685/suppl/GSE131685_RAW.tar"
        self.download_url_meta = None

        self.author = "Liao"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "kidney"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2020
        self.doi = "10.1038/s41597-019-0351-8"

        self.var_symbol_col = "names"
        self.var_ensembl_col = "ensembl"

    def _load(self):
        fn = os.path.join(self.data_dir, "GSE131685_RAW.tar")
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
