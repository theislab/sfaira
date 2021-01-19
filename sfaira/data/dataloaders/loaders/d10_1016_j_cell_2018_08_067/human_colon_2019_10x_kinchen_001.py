import anndata
import os
from typing import Union
import pandas as pd

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
        self.id = "human_colon_2019_10x_kinchen_001_10.1016/j.cell.2018.08.067"

        self.download = "https://data.humancellatlas.org/project-assets/project-matrices/f8aa201c-4ff1-45a4-890e-840d63459ca2.homo_sapiens.loom"
        self.download_meta = "private"

        self.author = "Simmons"
        self.doi = "10.1016/j.cell.2018.08.067"
        self.normalization = "raw"
        self.organ = "colon"  # ToDo: "lamina propria of mucosa of colon"
        self.organism = "human"
        self.protocol = "10x"
        self.year = 2019

        self.var_symbol_col = "names"
        self.var_ensembl_col = "Accession"

        self.obs_key_state_exact = "donor_organism.diseases.ontology_label"
        self.obs_key_healthy = self.obs_key_state_exact
        self.healthy_state_healthy = "normal"
        self.obs_key_cellontology_original = "celltype"

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

    def _load(self, fn=None):
        if fn is None:
            fn = [
                os.path.join(self.path, "human", "colon", "f8aa201c-4ff1-45a4-890e-840d63459ca2.homo_sapiens.loom"),
                os.path.join(self.path, "human", "colon", "uc_meta_data_stromal_with_donor.txt"),
                os.path.join(self.path, "human", "colon", "hc_meta_data_stromal_with_donor.txt")
            ]
        adata = anndata.read_loom(fn[0])
        ctuc = pd.read_csv(fn[1], sep="\t")
        cthealthy = pd.read_csv(fn[2], sep="\t")
        adata = adata[adata.obs["emptydrops_is_cell"] == "t"].copy()
        adata = adata[adata.X.sum(axis=1).flatten() >= 250].copy()
        uc = adata[adata.obs["donor_organism.diseases.ontology_label"] == "ulcerative colitis (disease)"].copy()
        bcuc = [i.split("-")[0] for i in ctuc["Barcode"]]
        seluc = []
        for i in uc.obs["barcode"]:
            seluc.append((uc.obs["barcode"].str.count(i).sum() == 1) and i in bcuc)
        uc = uc[seluc].copy()
        ctuc.index = [i.split("-")[0] for i in ctuc["Barcode"]]
        uc.obs["celltype"] = [ctuc.loc[i]["Cluster"] for i in uc.obs["barcode"]]
        uc.var = uc.var.reset_index().rename(columns={"index": "names"}).set_index("featurekey")
        healthy = adata[adata.obs["donor_organism.diseases.ontology_label"] == "normal"].copy()
        bchealthy = [i.split("-")[0] for i in cthealthy["Barcode"]]
        selhealthy = []
        for i in healthy.obs["barcode"]:
            selhealthy.append((healthy.obs["barcode"].str.count(i).sum() == 1) and i in bchealthy)
        healthy = healthy[selhealthy].copy()
        cthealthy.index = [i.split("-")[0] for i in cthealthy["Barcode"]]
        healthy.obs["celltype"] = [cthealthy.loc[i]["Cluster"] for i in healthy.obs["barcode"]]
        healthy.var = healthy.var.reset_index().rename(columns={"index": "names"}).set_index("featurekey")
        self.adata = healthy.concatenate(uc)
