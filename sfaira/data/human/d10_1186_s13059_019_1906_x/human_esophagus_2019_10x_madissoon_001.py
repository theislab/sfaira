import anndata
import os
from typing import Union
from .external import DatasetBase
import scipy.sparse


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
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_esophagus_2019_10x_madissoon_001_10.1186/s13059-019-1906-x"
        self.download = "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/oesophagus.cellxgene.h5ad"
        # Associated HCA project: https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7
        self.download_meta = None
        self.organ = "esophagus"
        self.sub_tissue = "esophagus"
        self.author = "Meyer"
        self.year = 2019
        self.doi = "10.1186/s13059-019-1906-x"
        self.protocol = "10x"
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.var_ensembl_col = 'gene_ids-HCATisStab7413619'
        self.obs_key_cellontology_original = 'Celltypes'

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

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "esophagus", "oesophagus.cellxgene.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                   .multiply(1/10000)
