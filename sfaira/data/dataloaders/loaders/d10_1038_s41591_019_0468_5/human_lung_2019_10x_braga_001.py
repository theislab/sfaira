import anndata
import os
from typing import Union
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_lung_2019_10x_braga_001_10.1038/s41591-019-0468-5"

        self.download = "https://covid19.cog.sanger.ac.uk/vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad"
        self.download_meta = None

        self.author = "Teichmann"
        self.doi = "10.1038/s41591-019-0468-5"
        self.healthy = True
        self.organ = "lung"  # ToDo: "alveoli, parenchyma"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2019
        self.normalization = "norm"

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "Ciliated 2": "Multiciliated lineage",
                "Luminal_Macrophages": "Macrophages",
                "Basal 1": "Basal",
                "Dendritic cells": "Dendritic cells",
                "Endothelial": "1_Endothelial",
                "Lymphatic": "Lymphatic EC",
                "Ciliated 1": "Multiciliated lineage",
                "Smooth muscle": "2_Smooth Muscle",
                "Type_1_alveolar": "AT1",
                "Neutrophils": "Monocytes",
                "Club": "Secretory",
                "Basal 2": "Basal",
                "B cells": "B cell lineage",
                "T and NK": "2_Lymphoid",
                "Mesothelium": "Mesothelium",
                "Mast cells": "Mast cells",
                "Fibroblasts": "2_Fibroblast lineage",
                "Type 2 alveolar": "AT2",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "lung", "vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
