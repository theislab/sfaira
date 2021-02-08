import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "lukassen20_lung_orig.processed.h5ad",
    "lukassen20_airway_orig.processed.h5ad"
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
        self.id = f"human_lung_2020_10x_lukassen_{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_" \
                  f"10.1101/2020.03.13.991455"

        self.download_url_data = f"https://covid19.cog.sanger.ac.uk/{self.sample_fn}"
        self.download_url_meta = None

        self.author = "Eils"
        self.doi = "10.1101/2020.03.13.991455"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        if self.sample_fn == "lukassen20_lung_orig.processed.h5ad":
            self.class_maps = {
                "0": {
                    "AT1": "AT1",
                    "AT2": "AT2",
                    "Ciliated": "Multiciliated lineage",
                    "Club": "Secretory",
                    "Endothelial": "1_Endothelial",
                    "Fibroblasts": "2_Fibroblast lineage",
                    "Immuno_TCells": "T cell lineage",
                    "Immuno_Monocytes": "Monocytes",
                    "LymphaticEndothelium": "Lymphatic EC",
                }
            }
        else:
            self.class_maps = {
                "0": {
                    "Basal_Mitotic": "Basal",
                    "Basal1": "Basal",
                    "Basal2": "Basal",
                    "Basal3": "Basal",
                    "Ciliated1": "Multiciliated lineage",
                    "Ciliated2": "Multiciliated lineage",
                    "Club": "Secretory",
                    "Fibroblast": "2_Fibroblast lineage",
                    "FOXN4": "Rare",
                    "Ionocyte": "Rare",
                    "Goblet": "Secretory",
                    "Secretory3": "Secretory",
                    "Secretory2": "Secretory",
                    "Secretory1": "Secretory",
                },
            }

    def _load(self, fn=None):
        base_path = os.path.join(self.path, "human", "lung")
        fn = os.path.join(base_path, self.sample_fn)
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["nCount_RNA"].values[:, None]))\
                                   .multiply(1 / 10000)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
