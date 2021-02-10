import anndata
import tarfile
import gzip
import scipy.io
import os
import pandas as pd
from typing import Union
from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "GSM3308545_NOD_08w_A",
    "GSM3308547_NOD_08w_C",
    "GSM3308548_NOD_14w_A",
    "GSM3308549_NOD_14w_B",
    "GSM3308550_NOD_14w_C",
    "GSM3308551_NOD_16w_A",
    "GSM3308552_NOD_16w_B",
    "GSM3308553_NOD_16w_C"
]


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def __init__(
            self,
            sample_fn: str,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = f"mouse_pancreas_2019_10x_thompson_{str(SAMPLE_FNS.index(sample_fn)).zfill(3)}_" \
                  f"10.1016/j.cmet.2019.01.021"

        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117770/suppl/GSE117770_RAW.tar"
        self.download_url_meta = f"private,{self.sample_fn}_annotation.csv"

        self.author = "Bhushan"
        self.doi = "10.1016/j.cmet.2019.01.021"
        self.healthy = False
        self.normalization = "raw"
        self.organ = "pancreas"
        self.organism = "mouse"
        self.protocol = "10X sequencing"
        self.state_exact = "diabetic"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "celltypes"

    def _load(self):
        with tarfile.open(os.path.join(self.doi_path, 'GSE117770_RAW.tar')) as tar:
            for member in tar.getmembers():
                if "_matrix.mtx.gz" in member.name and self.sample_fn in member.name:
                    name = "_".join(member.name.split("_")[:-1])
                    with gzip.open(tar.extractfile(member), "rb") as mm:
                        x = scipy.io.mmread(mm).T.tocsr()
                    obs = pd.read_csv(tar.extractfile(name + "_barcodes.tsv.gz"), compression="gzip", header=None,
                                      sep="\t", index_col=0)
                    obs.index.name = None
                    var = pd.read_csv(tar.extractfile(name + "_genes.tsv.gz"), compression="gzip", header=None,
                                      sep="\t")
                    var.columns = ["ensembl", "names"]
                    var.index = var["ensembl"].values
                    self.adata = anndata.AnnData(X=x, obs=obs, var=var)
        self.adata.var_names_make_unique()
        celltypes = pd.read_csv(os.path.join(self.doi_path, self.sample_fn + "_annotation.csv"), index_col=0)
        self.adata = self.adata[celltypes.index]
        self.adata.obs["celltypes"] = celltypes
