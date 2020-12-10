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
        self.id = "human_pancreas_2020_microwell_han_003_10.1038/s41586-020-2157-4"
        self.organ = 'Pancreas'
        self.sub_tissue = 'FetalPancreas'
        self.dev_stage = 'Fetus'
        self.download_website = 'https://figshare.com/articles/HCL_DGE_Data/7235471'
        self.download_website_meta = None
        self.annotated = True

        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'B cell': 'B cell',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'Basal cell': 'Basal cell',
                'CB CD34+': 'CB CD34+',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cell',
                'Endothelial cell (APC)': 'Endothelial cell',
                'Endothelial cell (endothelial to mesenchymal transition)': 'Endothelial cell',
                'Enterocyte progenitor': 'Enterocyte progenitor',
                'Erythroid cell': 'Erythroid cell',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Fetal Neuron': 'Neuron',
                'Fetal acinar cell': 'Acinar cell',
                'Fetal endocrine cell': 'Endocrine cell',
                'Fetal enterocyte ': 'Enterocyte',
                'Fetal epithelial progenitor': 'Epithelial progenitor',
                'Fetal fibroblast': 'Fibroblast',
                'Fetal mesenchymal progenitor': 'Mesenchymal Cell',
                'Fetal neuron': 'Neuron',
                'Fetal skeletal muscle cell': 'Skeletal muscle cell',
                'Fetal stromal cell': 'Stromal cell',
                'Fibroblast': 'Fibroblast',
                'Gastric endocrine cell': 'Gastric endocrine cell',
                'Immature sertoli cell (Pre-Sertoli cell)': 'Immature sertoli cell (Pre-Sertoli cell)',
                'Macrophage': 'Macrophage',
                'Mast cell': 'Mast cell',
                'Monocyte': 'Monocyte',
                'Neutrophil': 'Neutrophil',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Pancreas exocrine cell': 'Pancreas exocrine cell',
                'Primordial germ cell': 'Primordial germ cell',
                'Proliferating T cell': 'T cell',
                'Proximal tubule progenitor': 'Proximal tubule progenitor',
                'Sinusoidal endothelial cell': 'Endothelial cell',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "pancreas", "hcl_FetalPancreas_2.h5ad")
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
