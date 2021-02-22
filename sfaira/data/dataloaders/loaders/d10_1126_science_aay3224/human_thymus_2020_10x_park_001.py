import anndata
import os
from typing import Union
import numpy as np

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
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/park20.processed.h5ad"
        self.download_url_meta = None

        self.author = "Park"
        self.doi = "10.1126/science.aay3224"
        self.healthy = True
        self.normalization = "norm"
        self.organ = "thymus"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "Anno_level_fig1"

        self.class_maps = {
            "0": {
                "B_memory": "B_memory",
                "B_naive": "B_naive",
                "B_plasma": "B_plasma",
                "B_pro/pre": "B_pro/pre",
                "CD4+T": "CD4+T",
                "CD4+Tmem": "CD4+Tmem",
                "CD8+T": "CD8+T",
                "CD8+Tmem": "CD8+Tmem",
                "CD8αα": "CD8αα",
                "DC1": "DC1",
                "DC2": "DC2",
                "DN": "DN",
                "DP": "DP",
                "ETP": "ETP",
                "Endo": "Endo",
                "Epi_GCM2": "Epi_GCM2",
                "Ery": "Ery",
                "Fb_1": "Fb_1",
                "Fb_2": "Fb_2",
                "Fb_cycling": "Fb_cycling",
                "ILC3": "ILC3",
                "Lymph": "Lymph",
                "Mac": "Mac",
                "Mast": "Mast",
                "Mgk": "Mgk",
                "Mono": "Mono",
                "NK": "NK",
                "NKT": "NKT",
                "NMP": "NMP",
                "T(agonist)": "T(agonist)",
                "TEC(myo)": "TEC(myo)",
                "TEC(neuro)": "TEC(neuro)",
                "Treg": "Treg",
                "VSMC": "VSMC",
                "aDC": "aDC",
                "cTEC": "cTEC",
                "mTEC(I)": "mTEC(I)",
                "mTEC(II)": "mTEC(II)",
                "mTEC(III)": "mTEC(III)",
                "mTEC(IV)": "mTEC(IV)",
                "mcTEC": "mcTEC",
                "pDC": "pDC",
                "αβT(entry)": "alpha_beta_T(entry)",
                "γδT": "gamma_delta_T",
            },
        }

    def _load(self):
        fn = os.path.join(self.data_dir, "park20.processed.h5ad")
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)

        return adata
