import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class. This dataloader only provides the subset of the published sata which has been made available through the
    covid-19 Cell Atlas.

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
        self.id = "human_colon_2019_10x_smilie_001_10.1016/j.cell.2019.06.029"

        self.download = "https://covid19.cog.sanger.ac.uk/smillie19_epi.processed.h5ad"
        self.download_meta = None

        self.author = "Regev"
        self.doi = "10.1016/j.cell.2019.06.029"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "colon"  # ToDo: "colonic epithelium"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "Cycling TA": "Cycling TA",
                "TA 1": "TA 1",
                "TA 2": "TA 2",
                "Immature Enterocytes 2": "Immature Enterocytes 2",
                "Immature Enterocytes 1": "Immature Enterocytes 1",
                "Enterocyte Progenitors": "Enterocyte Progenitors",
                "Immature Goblet": "Immature Goblet",
                "Enterocytes": "Enterocytes",
                "Secretory TA": "Secretory TA",
                "Best4+ Enterocytes": "Best4+ Enterocytes",
                "CD8+ IELs": "CD8+ IELs",
                "Goblet": "Goblet cells",
                "Stem": "Stem cells",
                "Tuft": "Tuft",
                "Follicular": "Follicular",
                "Enteroendocrine": "Enteroendocrine cells",
                "Plasma": "Plasma Cells",
                "CD4+ Memory": "CD4+ Memory",
                "CD8+ LP": "CD8+ LP",
                "CD69- Mast": "CD69- Mast",
                "Macrophages": "Macrophage",
                "GC": "Glial cells",
                "Cycling B": "B cell cycling",
                "CD4+ Activated Fos-hi": "CD4+ T Activated Fos-hi",
                "CD4+ Activated Fos-lo": "CD4+ T Activated Fos-lo",
                "NKs": "NK",
                "Cycling T": "Cycling T",
                "M cells": "M cells",
                "CD69+ Mast": "CD69+ Mast",
                "MT-hi": "MT-hi",
                "CD8+ IL17+": "CD8+ IL17+",
                "CD4+ PD1+": "CD4+ PD1+",
                "DC2": "DC2",
                "Treg": "Treg",
                "ILCs": "ILC",
                "DC1": "DC1",
                "WNT2B+ Fos-lo 1": "WNT2B+ Fos-lo 1",
                "WNT5B+ 2": "WNT5B+ 2",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "colon", "smillie19_epi.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["n_counts"].values[:, None]))\
                                   .multiply(1 / 10000)
