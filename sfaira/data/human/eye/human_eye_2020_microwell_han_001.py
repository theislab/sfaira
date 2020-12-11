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
        self.id = "human_eye_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'Eye'
        self.sub_tissue = 'FetalEyes'
        self.dev_stage = 'Fetus'
        self.download = 'https://figshare.com/articles/HCL_DGE_Data/7235471'
        self.download_meta = None
        self.annotated = True

        self.class_maps = {
            "0": {
                'Fetal neuron': 'Fetal neuron',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Erythroid cell': 'Erythroid cell',
                'Primordial germ cell': 'Primordial germ cell',
                'Endothelial cell': 'Endothelial cell',
                'Fetal skeletal muscle cell': 'Fetal skeletal muscle cell',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Fetal fibroblast': 'Fibroblast',
                'Fetal Neuron': 'Fetal neuron',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'Dendritic cell': 'Dendritic cell',
                'Fetal endocrine cell': 'Fetal endocrine cell',
                'Macrophage': 'Macrophage',
                'T cell': 'T cell',
                'Basal cell': 'Basal cell',
                'Gastric endocrine cell': 'Gastric endocrine cell',
                'Goblet cell': 'Goblet cell',
                'Epithelial cell (intermediated)': 'Epithelial cell (intermediated)',
                'Stratified epithelial cell': 'Stratified epithelial cell',
                'CB CD34+': 'CB CD34_pos',
                'hESC': 'hESC'
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "eye", "hcl_FetalEyes_1.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = 'Guo'
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = 2020
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = '10.1038/s41586-020-2157-4'
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = "microwell"
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = self.species
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.download_meta] = self.download_meta
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'
        self.adata.uns[self._ADATA_IDS_SFAIRA.dev_stage] = self.dev_stage

        self._convert_and_set_var_names(symbol_col="names", ensembl_col="ensembl")

