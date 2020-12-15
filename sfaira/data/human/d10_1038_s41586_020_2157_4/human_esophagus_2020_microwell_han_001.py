import os
from typing import Union
from .external import DatasetHcl


class Dataset(DatasetHcl):
    """
    This is a dataloader for a the Human Cell Landscape dataset (Han et al. 2020. doi: 10.1038/s41586-020-2157-4).

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
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)        self.id = "human_esophagus_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'Esophagus'
        self.sub_tissue = 'AdultEsophagus'
        self.dev_stage = 'Adult'
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
        self._load_hcl(fn=fn, sample_id="AdultEsophagus_1")
