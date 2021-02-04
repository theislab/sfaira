import anndata
import numpy as np
import os
import pandas
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

        self.download = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117770"
        self.download_meta = "private"

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

    def _load(self, fn=None):
        path_base = os.path.join(self.path, "mouse", "pancreas")
        celltypes = pandas.read_csv(os.path.join(path_base, self.sample_fn + "_annotation.csv"), index_col=0)

        self.adata = anndata.read_mtx(os.path.join(path_base, self.sample_fn + "_matrix.mtx.gz")).transpose()
        self.adata.var_names = np.genfromtxt(os.path.join(path_base, self.sample_fn + "_genes.tsv.gz"), dtype=str)[:, 1]
        self.adata.obs_names = np.genfromtxt(os.path.join(path_base, self.sample_fn + "_barcodes.tsv.gz"), dtype=str)
        self.adata.var_names_make_unique()
        self.adata = self.adata[celltypes.index]
        self.adata.obs["celltypes"] = celltypes
