import anndata
import os
from typing import Union
from .external import DatasetBase


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
        self.id = "human_esophagus_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'Esophagus'
        self.sub_tissue = 'AdultEsophagus'
        self.dev_stage = 'Adult'
        self.download = 'https://figshare.com/articles/HCL_DGE_Data/7235471'
        self.download_meta = None
        self.annotated = True

        self.class_maps = {
            "0": {
                'Fibroblast': 'Fibroblast',
                'Basal cell': 'Basal cell',
                'Stratified epithelial cell': 'Stratified epithelial cell',
                'Endothelial cell (APC)': 'Endothelial cell (APC)',
                'Macrophage': 'Macrophage',
                'B cell': 'B cell',
                'T cell': 'T cell',
                'Dendritic cell': 'Dendritic cell',
                'Mast cell': 'Mast cell',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'Stromal cell': 'Stromal cell',
                'Monocyte': 'Monocyte',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Endothelial cell': 'Endothelial cell',
                'Neutrophil': 'Neutrophil',
                'Endothelial cell (endothelial to mesenchymal transition)': 'Endothelial cell (endothelial to mesenchymal transition)',
                'Fetal stromal cell': 'Fetal stromal cell',
                'CB CD34+': 'CB CD34+',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Gastric endocrine cell': 'Gastric endocrine cell',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Sinusoidal endothelial cell': 'Sinusoidal endothelial cell',
                'Loop of Henle': 'Loop of Henle',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "esophagus", "hcl_AdultEsophagus_1.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = 'Guo'
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = 2020
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = '10.1038/s41586-020-2157-4'
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = "microwell"
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'
        self.adata.uns[self._ADATA_IDS_SFAIRA.dev_stage] = self.dev_stage

        self._convert_and_set_var_names(symbol_col="names", ensembl_col="ensembl")
