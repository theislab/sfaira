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
        self.id = "human_lung_2019_10x_braga_002_10.1038/s41591-019-0468-5"

        self.download = "https://covid19.cog.sanger.ac.uk/vieira19_Bronchi_anonymised.processed.h5ad"
        self.download_meta = None

        self.author = "Teichmann"
        self.doi = "10.1038/s41591-019-0468-5"
        self.healthy = True
        self.normalization = "norm"
        self.organ = "lung"  # ToDo "bronchi"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "Ciliated 1": "Multiciliated lineage",
                "Club": "Secretory",
                "Ciliated 2": "Multiciliated lineage",
                "Ionocytes": "Rare",
                "Basal 2": "Basal",
                "Goblet_1": "Secretory",
                "Goblet 2": "Secretory",
                "Basal 1": "Basal",
                "Dendritic cells": "Dendritic cells",
                "B cells": "B cell lineage",
                "Luminal_Macrophages": "Macrophages",
                "Neutrophils": "Monocytes",
                "Endothelial": "1_Endothelial",
                "Smooth muscle": "2_Smooth Muscle",
                "T and NK": "2_Lymphoid",
                "Fibroblasts": "2_Fibroblast lineage",
                "Lymphatic": "Lymphatic EC",
                "Mast cells": "Mast cells",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "lung", "vieira19_Bronchi_anonymised.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
