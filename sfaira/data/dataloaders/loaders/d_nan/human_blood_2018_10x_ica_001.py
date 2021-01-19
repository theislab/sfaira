import anndata
import os
from typing import Union
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):

        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_blood_2018_10x_ica_001_unknown"

        self.download_url_data = "https://data.humancellatlas.org/project-assets/project-matrices/cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom"
        self.download_url_meta = None

        self.author = "Regev"
        self.doi = "d_nan"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "blood"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "index"
        self.var_ensembl_col = "Accession"

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "blood", "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom")
        self.adata = anndata.read_loom(fn)
        idx = np.logical_and((self.adata.obs["derived_organ_parts_label"] == "umbilical cord blood").values,
                             (self.adata.obs["emptydrops_is_cell"] == "t").values)
        self.adata = self.adata[idx].copy()
