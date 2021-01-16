import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

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
        self.id = "human_lung_2020_10x_lukassen_001_10.1101/2020.03.13.991455"

        self.download = "https://covid19.cog.sanger.ac.uk/lukassen20_lung_orig.processed.h5ad"
        self.download_meta = None

        self.author = "Eils"
        self.doi = "10.1101/2020.03.13.991455"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "Ciliated": "Multiciliated lineage",
                "Endothelial": "1_Endothelial",
                "AT2": "AT2",
                "LymphaticEndothelium": "Lymphatic EC",
                "Fibroblasts": "2_Fibroblast lineage",
                "Club": "Secretory",
                "Immuno_TCells": "T cell lineage",
                "Immuno_Monocytes": "Monocytes",
                "AT1": "AT1"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "lung", "lukassen20_lung_orig.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["nCount_RNA"].values[:, None]))\
                                   .multiply(1 / 10000)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
