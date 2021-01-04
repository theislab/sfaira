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
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_pancreas_2020_microwell_han_003_10.1038/s41586-020-2157-4"
        self.organ = 'Pancreas'
        self.sub_tissue = 'FetalPancreas'
        self.dev_stage = 'Fetus'
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
        self._load_hcl(fn=fn, sample_id="FetalPancreas_2")
