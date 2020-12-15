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
        self.id = "human_brain_2020_microwell_han_006_10.1038/s41586-020-2157-4"
        self.organ = 'brain'
        self.sub_tissue = 'AdultCerebellum'
        self.dev_stage = 'Adult'
        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'Astrocyte': 'Astrocyte',
                'B cell': 'B cell',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'CB CD34+': 'CB CD34+',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cells',
                'Endothelial cell (APC)': 'Endothelial cells',
                'Erythroid cell': 'Erythroid cell',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Fetal Neuron': 'Fetal Neuron',
                'Fetal endocrine cell': 'Fetal endocrine cell',
                'Fetal enterocyte ': 'Fetal enterocyte ',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Fetal neuron': 'Fetal Neuron',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Fibroblast': 'Fibroblast',
                'Gastric endocrine cell': 'Gastric endocrine cell',
                'Goblet cell': 'Goblet cell',
                'Macrophage': 'Macrophage',
                'Monocyte': 'Monocyte',
                'Neutrophil': 'Neutrophil',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Oligodendrocyte': 'Oligodendrocytes',
                'Primordial germ cell': 'Primordial germ cell',
                'Sinusoidal endothelial cell': 'Sinusoidal endothelial cell',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
                'hESC': 'Neuronal stem cells'
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="AdultCerebellum_1")
