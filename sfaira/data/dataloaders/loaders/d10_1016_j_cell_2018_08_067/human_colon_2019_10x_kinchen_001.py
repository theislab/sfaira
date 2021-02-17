import os
from typing import Union
import pandas as pd
import anndata as ad
import scipy.sparse
import numpy as np

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "HC",
    "UC",
]


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def __init__(
            self,
            sample_fn: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = f"human_colon_2019_10x_kinchen_{str(SAMPLE_FNS.index(sample_fn)+1).zfill(3)}_10.1016/j.cell.2018.08.067"

        self.download_url_data = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114374&format=file&" \
                                 f"file=GSE114374%5FHuman%5F{sample_fn}%5Fexpression%5Fmatrix%2Etxt%2Egz"
        self.download_url_meta = f"private,{sample_fn.lower()}_meta_data_stromal_with_donor.txt"

        self.author = "Simmons"
        self.doi = "10.1016/j.cell.2018.08.067"
        self.normalization = "norm"
        self.organ = "lamina propria of mucosa of colon"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.year = 2019

        self.var_symbol_col = "index"
        self.obs_key_state_exact = "state_exact"
        self.obs_key_healthy = self.obs_key_state_exact
        self.healthy_state_healthy = "healthy colon"
        self.obs_key_cellontology_original = "Cluster"
        self.obs_key_age = "Age"
        self.obs_key_sex = "Sex"

        self.class_maps = {
            "0": {
                "Endothelial 1": "Endothelial",
                "Endothelial 2": "Endothelial",
                "Glial": "Glial cells",
                "Myofibroblasts": "Myofibroblasts",
                "Pericyte 1": "Pericytes",
                "Pericyte 2": "Pericytes",
                "Pericytes": "Pericytes",
                "Plasma Cells": "Plasma Cells",
                "Smooth Muscle": "Smooth Muscle",
                "Stromal 1": "Stromal",
                "Stromal 2a": "Stromal",
                "Stromal 2b": "Stromal",
                "Stromal 3": "Stromal",
                "Stromal 4": "Stromal",
            },
        }

    def _load(self):
        fn = [
            os.path.join(self.data_dir, f"GSE114374_Human_{self.sample_fn}_expression_matrix.txt.gz"),
            os.path.join(self.data_dir, f"{self.sample_fn.lower()}_meta_data_stromal_with_donor.txt"),
        ]
        matrix = pd.read_csv(fn[0], sep="\t")
        obs = pd.read_csv(fn[1], sep="\t", index_col=3)
        self.adata = ad.AnnData(matrix.T)
        self.adata.X = scipy.sparse.csc_matrix(np.expm1(self.adata.X))
        self.adata.obs = obs
        self.adata.obs['state_exact'] = "healthy colon" if self.sample_fn == "HC" else "ulcerative colitis"
        s_dict = {"F": "female", "M": "male"}
        self.adata.obs['Sex'] = [s_dict[i] for i in self.adata.obs['Sex']]
