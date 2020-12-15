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
        super().__init__(path=path, meta_path=meta_path, **kwargs)
        self.id = "human_spleen_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'Spleen'
        self.sub_tissue = 'AdultSpleen'
        self.dev_stage = 'Adult'
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
        self._load_hcl(fn=fn, sample_id="AdultSpleenParenchyma_1")
