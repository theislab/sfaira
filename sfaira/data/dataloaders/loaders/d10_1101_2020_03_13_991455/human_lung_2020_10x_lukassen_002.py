import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

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
        self.id = "human_lung_2020_10x_lukassen_002_10.1101/2020.03.13.991455"

        self.download = "https://covid19.cog.sanger.ac.uk/lukassen20_airway_orig.processed.h5ad"
        self.download_meta = None

        self.author = "Eils"
        self.doi = "10.1101/2020.03.13.991455"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"  # ToDo: "bronchial epithelial cells"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "Secretory3": "Secretory",
                "Ciliated1": "Multiciliated lineage",
                "Goblet": "Secretory",
                "Ciliated2": "Multiciliated lineage",
                "Club": "Secretory",
                "Secretory2": "Secretory",
                "FOXN4": "Rare",
                "Basal1": "Basal",
                "Secretory1": "Secretory",
                "Fibroblast": "2_Fibroblast lineage",
                "Ionocyte": "Rare",
                "Basal3": "Basal",
                "Basal_Mitotic": "Basal",
                "Basal2": "Basal",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "lung", "lukassen20_airway_orig.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["nCount_RNA"].values[:, None]))\
                                   .multiply(1 / 10000)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
