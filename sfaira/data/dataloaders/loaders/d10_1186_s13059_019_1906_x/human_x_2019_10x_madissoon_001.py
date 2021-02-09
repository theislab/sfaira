import anndata
import os
from typing import Union
import scipy.sparse

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "madissoon19_lung.processed.h5ad",
    "oesophagus.cellxgene.h5ad",
    "spleen.cellxgene.h5ad",
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
        organ = "lung parenchyma" if self.sample_fn == "madissoon19_lung.processed.h5ad" else \
            "esophagus" if self.sample_fn == "oesophagus.cellxgene.h5ad" else "spleen"
        self.id = f"human_{"".join(organ.split(" "))}_2019_10x_madissoon_{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_" \
                  f"10.1186/s13059-019-1906-x"

        if self.sample_fn == "madissoon19_lung.processed.h5ad":
            "https://covid19.cog.sanger.ac.uk/madissoon19_lung.processed.h5ad"
            self.var_ensembl_col = "gene.ids.HCATisStab7509734"
        elif self.sample_fn == "oesophagus.cellxgene.h5ad":
            self.download = "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/oesophagus.cellxgene.h5ad"
            # Associated HCA project: https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7
            self.var_ensembl_col = "gene_ids-HCATisStab7413619"
        else:
            "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/spleen.cellxgene.h5ad"
            self.var_ensembl_col = "gene_ids-HCATisStab7463846"

        self.author = "Meyer"
        self.doi = "10.1186/s13059-019-1906-x"
        self.healthy = True
        self.normalization = "raw"
        self.organ = organ
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "Celltypes"

        if self.sample_fn == "madissoon19_lung.processed.h5ad":
            self.class_maps = {
                "0": {
                    "T_CD4": "T cell lineage",
                    "Mast_cells": "Mast cells",
                    "Monocyte": "Monocytes",
                    "Blood_vessel": "2_Blood vessels",
                    "Ciliated": "Multiciliated lineage",
                    "Macrophage_MARCOneg": "Macrophages",
                    "DC_plasmacytoid": "Dendritic cells",
                    "DC_1": "Dendritic cells",
                    "Muscle_cells": "2_Smooth Muscle",
                    "Macrophage_MARCOpos": "Macrophages",
                    "T_cells_Dividing": "T cell lineage",
                    "DC_Monocyte_Dividing": "Dendritic cells",
                    "B_cells": "B cell lineage",
                    "T_CD8_CytT": "T cell lineage",
                    "NK_Dividing": "Innate lymphoid cells",
                    "T_regulatory": "T cell lineage",
                    "DC_2": "Dendritic cells",
                    "Alveolar_Type2": "AT2",
                    "Plasma_cells": "B cell lineage",
                    "NK": "Innate lymphoid cells",
                    "Alveolar_Type1": "AT1",
                    "Fibroblast": "2_Fibroblast lineage",
                    "DC_activated": "Dendritic cells",
                    "Macrophage_Dividing": "Macrophages",
                    "Lymph_vessel": "Lymphatic EC",
                },
            }
        elif self.sample_fn == "oesophagus.cellxgene.h5ad":
            self.class_maps = {
                "0": {
                    "B_CD27neg": "B_CD27neg",
                    "B_CD27pos": "B_CD27pos",
                    "Blood_vessel": "Blood_vessel",
                    "Dendritic_Cells": "Dendritic cell",
                    "Epi_basal": "Basal cell",
                    "Epi_dividing": "Epi_dividing",
                    "Epi_stratified": "Stratified epithelial cell",
                    "Epi_suprabasal": "Epi_suprabasal",
                    "Epi_upper": "Epi_upper",
                    "Glands_duct": "Glands_duct",
                    "Glands_mucous": "Glands_mucous",
                    "Lymph_vessel": "Lymph_vessel",
                    "Mast_cell": "Mast cell",
                    "Mono_macro": "Mono_macro",
                    "NK_T_CD8_Cytotoxic": "NK_T_CD8_Cytotoxic",
                    "Stroma": "Stromal cell",
                    "T_CD4": "T_CD4",
                    "T_CD8": "T_CD8",
                },
            }
        else:
            self.class_maps = {
                "0": {
                    "B_Hypermutation": "B_Hypermutation",
                    "B_T_doublet": "B_T_doublet",
                    "B_follicular": "B_follicular",
                    "B_mantle": "B_mantle",
                    "CD34_progenitor": "CD34_progenitor",
                    "DC_1": "DC_1",
                    "DC_2": "DC_2",
                    "DC_activated": "DC_activated",
                    "DC_plasmacytoid": "DC_plasmacytoid",
                    "ILC": "ILC",
                    "Macrophage": "Macrophage",
                    "Monocyte": "Monocyte",
                    "NK_CD160pos": "NK_CD160pos",
                    "NK_FCGR3Apos": "NK_FCGR3Apos",
                    "NK_dividing": "NK_dividing",
                    "Plasma_IgG": "Plasma_IgG",
                    "Plasma_IgM": "Plasma_IgM",
                    "Plasmablast": "Plasmablast",
                    "Platelet": "Platelet",
                    "T_CD4_conv": "T_CD4_conv",
                    "T_CD4_fh": "T_CD4_fh",
                    "T_CD4_naive": "T_CD4_naive",
                    "T_CD4_reg": "T_CD4_reg",
                    "T_CD8_CTL": "T_CD8_CTL",
                    "T_CD8_MAIT": "T_CD8_MAIT",
                    "T_CD8_activated": "T_CD8_activated",
                    "T_CD8_gd": "T_CD8_gd",
                    "T_cell_dividing": "Proliferating T cell",
                },
            }

    def _load(self, fn=None):
        base_path = os.path.join(self.path, "human", self.organ)
        fn = os.path.join(base_path, self.sample_fn)
        self.adata = anndata.read(fn)
        if self.sample_fn == "oesophagus.cellxgene.h5ad" or self.sample_fn == "spleen.cellxgene.h5ad":
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["n_counts"].values[:, None]))\
                                       .multiply(1 / 10000)
        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
