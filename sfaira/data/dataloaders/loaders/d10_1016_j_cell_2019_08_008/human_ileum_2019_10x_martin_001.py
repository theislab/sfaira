import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

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
        self.id = "human_ileum_2019_10x_martin_001_10.1016/j.cell.2019.08.008"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/martin19.processed.h5ad"
        self.download_url_meta = None

        self.author = "Kenigsberg"
        self.doi = "10.1016/j.cell.2019.08.008"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "ileum"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019
        self.var_symbol_col = "index"
        self.var_ensembl_col = "gene_ids"
        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "T cells": "T cells",
                "Plasma cells": "Plasma Cells",
                "B cells": "B cells",
                "MNP": "MNP",
                "ILC": "ILC",
                "Enterocytes": "Enterocytes",
                "Fibs": "Fibroblasts",
                "CD36+ endothelium": "CD36+ endothelium",
                "Progenitors": "Progenitors",
                "Goblets": "Goblet cells",
                "Glial cells": "Glial cells",
                "Cycling": "Cycling",
                "ACKR1+ endothelium": "ACKR1+ endothelium",
                "Pericytes": "Pericytes",
                "Lymphatics": "Lymphatics",
                "Mast cells": "Mast cells",
                "SM": "Smooth muscle cell",
                "TA": "TA",
                "Paneth cells": "Paneth cells",
                "Enteroendocrines": "Enteroendocrine cells",
            },
        }

    def _load(self):
        fn = os.path.join(self.data_dir, "martin19.processed.h5ad")
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None]))\
                                   .multiply(1 / 10000)
        adata = adata[adata.obs["CellType"] != "Doublets"].copy()

        return adata
