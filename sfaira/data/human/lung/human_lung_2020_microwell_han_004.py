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
        self.id = "human_lung_2020_microwell_han_004_10.1038/s41586-020-2157-4"
        self.organ = 'lung'
        self.sub_tissue = 'AdultLung'
        self.dev_stage = 'Adult'
        self.download = 'https://figshare.com/articles/HCL_DGE_Data/7235471'
        self.download_meta = None
        self.annotated = True

        self.class_maps = {
            "0": {
                'AT2 cell': 'AT2',
                'Antigen presenting cell (RPS high)': 'unknown',
                'B cell': 'B cell lineage',
                'B cell (Plasmocyte)': 'B cell lineage',
                'Basal cell': 'Basal',
                'CB CD34+': 'Fetal airway progenitors',
                'Chondrocyte': '1_Stroma',
                'Dendritic cell': 'Dendritic cells',
                'Endothelial cell': '1_Endothelial',
                'Endothelial cell (APC)': '1_Endothelial',
                'Endothelial cell (endothelial to mesenchymal transition)': '1_Endothelial',
                'Enterocyte progenitor': '1_Epithelial',
                'Epithelial cell': '1_Epithelial',
                'Epithelial cell (intermediated)': '1_Epithelial',
                'Erythroid cell': 'Erythrocytes',
                'Erythroid progenitor cell (RP high)': 'Erythrocytes',
                'Fasciculata cell': 'unknown',
                'Fetal Neuron': 'unknown',
                'Fetal chondrocyte': '1_Stroma',
                'Fetal endocrine cell': 'unknown',
                'Fetal enterocyte ': '1_Epithelial',
                'Fetal epithelial progenitor': '1_Epithelial',
                'Fetal fibroblast': 'Fibroblasts',
                'Fetal mesenchymal progenitor': '1_Stroma',
                'Fetal neuron': 'unknown',
                'Fetal skeletal muscle cell': 'unknown',
                'Fetal stromal cell': '1_Stroma',
                'Fibroblast': 'Fibroblasts',
                'Gastric endocrine cell': 'unknown',
                'Goblet cell': 'Secretory',
                'Kidney intercalated cell': 'unknown',
                'Loop of Henle': 'unknown',
                'M2 Macrophage': 'Macrophages',
                'Macrophage': 'Macrophages',
                'Mast cell': 'Mast cells',
                'Mesothelial cell': 'Mast cells',
                'Monocyte': 'Monocytes',
                'Myeloid cell': '2_Myeloid',
                'Neutrophil': 'Neutrophilic',
                'Neutrophil (RPS high)': 'Neutrophilic',
                'Primordial germ cell': 'unknown',
                'Proliferating T cell': 'T cell lineage',
                'Proximal tubule progenitor': 'unknown',
                'Sinusoidal endothelial cell': '1_Endothelial',
                'Smooth muscle cell': '2_Smooth Muscle',
                'Stratified epithelial cell': '1_Epithelial',
                'Stromal cell': '1_Stroma',
                'T cell': 'T cell lineage',
                'Ventricle cardiomyocyte': '1_Stroma',
                'hESC': 'Fetal airway progenitors',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "lung", "hcl_AdultLung_1.h5ad")
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

