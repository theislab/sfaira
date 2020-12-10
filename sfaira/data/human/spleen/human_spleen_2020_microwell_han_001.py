import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata


class Dataset(DatasetBase):
    """
    This is a dataloader for a the Human Cell Landscape dataset (Han et al. 2020. doi: 10.1038/s41586-020-2157-4).
    In order to obtain the required preprocessed datafiles, please use the notebook provided in this repository under:
    sfaira/data/download_scripts/get_and_preprocess_HumanCellLandscape.ipynb

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
        self.id = "human_spleen_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'Spleen'
        self.sub_tissue = 'AdultSpleen'
        self.dev_stage = 'Adult'
        self.download_website = 'https://figshare.com/articles/HCL_DGE_Data/7235471'
        self.download_website_meta = None
        self.annotated = True

        self.class_maps = {
            "0": {
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'Neutrophil': 'Neutrophil',
                'Endothelial cell (APC)': 'Endothelial cell (APC)',
                'B cell': 'B cell',
                'Macrophage': 'Macrophage',
                'T cell': 'T cell',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Dendritic cell': 'Dendritic cell',
                'CB CD34+': 'CB CD34+',
                'Erythroid cell': 'Erythroid cell',
                'Monocyte': 'Monocyte',
                'Endothelial cell': 'Endothelial cell',
                'Sinusoidal endothelial cell': 'Sinusoidal endothelial cell',
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Proliferating T cell': 'Proliferating T cell',
                'Fibroblast': 'Fibroblast',
                'Stromal cell': 'Stromal cell',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Mast cell': 'Mast cell',
                'Smooth muscle cell': 'Smooth muscle cell',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "spleen", "hcl_AdultSpleenParenchyma_1.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Guo'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2020
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = '10.1038/s41586-020-2157-4'
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = "microwell"
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'
        self.adata.uns["dev_stage"] = self.dev_stage

        self._convert_and_set_var_names(symbol_col=ADATA_IDS_SFAIRA.gene_id_names, ensembl_col=ADATA_IDS_SFAIRA.gene_id_ensembl)
