import anndata
import tarfile
import gzip
import scipy.io
import os
import pandas as pd
from typing import Union
from sfaira.data import DatasetBase


class Dataset_d10_1016_j_cmet_2019_01_021(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117770/suppl/GSE117770_RAW.tar"
        self.download_url_meta = "private"

        self.author = "Bhushan"
        self.doi = "10.1016/j.cmet.2019.01.021"
        self.healthy = False
        self.normalization = "raw"
        self.organ = "pancreas"
        self.organism = "mouse"
        self.protocol = "10x"
        self.state_exact = "diabetic"
        self.year = 2019

        self.var_symbol_col = "index"

        self.class_maps = {
            "0": {
                "acinar": "pancreatic acinar cell",
                "ductal": "pancreatic ductal cell",
                "leukocyte": "leukocyte",
                "T cell(Pancreas)": "t cell",
                "B cell(Pancreas)": "b cell",
                "beta": "pancreatic B cell",
                "alpha": "pancreatic A cell",
                "delta": "pancreatic D cell",
                "pp": "pancreatic PP cell",
                "smooth_muscle": "smooth muscle cell",
                "stellate cell": "pancreatic stellate cell",
                "fibroblast": "stromal cell",
                "endothelial": "endothelial cell"
            },
        }

    def _load_generalized(self, data_path, sample_name, path_meta):

        with tarfile.open(os.path.join(data_path, 'GSE117770_RAW.tar')) as tar:
            for member in tar.getmembers():
                if "_matrix.mtx.gz" in member.name and sample_name in member.name:
                    name = "_".join(member.name.split("_")[:-1])
                    with gzip.open(tar.extractfile(member), "rb") as mm:
                        X = scipy.io.mmread(mm).T.tocsr()
                    obs = pd.read_csv(tar.extractfile(name + "_barcodes.tsv.gz"), compression="gzip", header=None,
                                      sep="\t", index_col=0)
                    obs.index.name = None
                    var = pd.read_csv(tar.extractfile(name + "_genes.tsv.gz"), compression="gzip", header=None,
                                      sep="\t")
                    var.columns = ["ensembl", "names"]
                    var.index = var["ensembl"].values
                    self.adata = anndata.AnnData(X=X, obs=obs, var=var)
        self.adata.var_names_make_unique()
        celltypes = pd.read_csv(path_meta, index_col=0)
        self.adata = self.adata[celltypes.index]
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = celltypes
