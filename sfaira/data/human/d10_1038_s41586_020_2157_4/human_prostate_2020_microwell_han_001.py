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
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.id = "human_prostate_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'prostate'
        self.sub_tissue = 'AdultProstate'
        self.dev_stage = 'Adult'
        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'Basal cell': 'Basal cell',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cell',
                'Endothelial cell (APC)': 'Endothelial cell',
                'Endothelial cell (endothelial to mesenchymal transition)': 'Endothelial cell',
                'Enterocyte progenitor': 'Enterocyte progenitor',
                'Epithelial cell (intermediated)': 'Epithelial cell (intermediated)',
                'Fasciculata cell': 'Fasciculata cell',
                'Fetal enterocyte': 'Fetal enterocyte',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Gastric endocrine cell': 'Gastric endocrine cell',
                'Goblet cell': 'Goblet cell',
                'Macrophage': 'Macrophage',
                'Monocyte': 'Monocyte',
                'Primordial germ cell': 'Primordial germ cell',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stratified epithelial cell': 'Stratified epithelial cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="AdultProstate_1")
