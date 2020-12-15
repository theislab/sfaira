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
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)        self.id = "human_skin_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'skin'
        self.sub_tissue = 'FetalSkin'
        self.dev_stage = 'Fetus'
        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'B cell': 'B cell',
                'Basal cell': 'Basal cell',
                'CB CD34+': 'CB CD34+',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cell',
                'Endothelial cell (APC)': 'Endothelial cell (APC)',
                'Epithelial cell': 'Epithelial cell',
                'Erythroid cell': 'Erythroid cell',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Fetal Neuron': 'Fetal Neuron',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Fetal fibroblast': 'Fetal fibroblast',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Fetal skeletal muscle cell': 'Fetal skeletal muscle cell',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Fibroblast': 'Fibroblast',
                'Kidney intercalated cell': 'Kidney intercalated cell',
                'Macrophage': 'Macrophage',
                'Mast cell': 'Mast cell',
                'Monocyte': 'Monocyte',
                'Neutrophil': 'Neutrophil',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Primordial germ cell': 'Primordial germ cell',
                'Proliferating T cell': 'Proliferating T cell',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
                'hESC': 'hESC',
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="FetalSkin_2")
