import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata
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
        self.id = "human_esophagus_2019_10x_madissoon_001_10.1101/741405"
        self.download_website = "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/oesophagus.cellxgene.h5ad"
        # Associated HCA project: https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7
        self.download_website_meta = None
        self.organ = "esophagus"
        self.sub_tissue = "esophagus"
        self.has_celltypes = True

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
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "esophagus", "oesophagus.cellxgene.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = "Meyer"
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1101/741405"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = "10x"
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['Celltypes']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col='gene_ids-HCATisStab7413619',
                                        new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)