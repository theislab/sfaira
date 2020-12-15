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
        self.id = "human_colon_2020_microwell_han_003_10.1038/s41586-020-2157-4"
        self.organ = 'Colon'
        self.sub_tissue = 'AdultTransverseColon'
        self.dev_stage = 'Adult'
        self.class_maps = {
            "0": {
                'Enterocyte progenitor': 'Enterocyte Progenitors',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'Enterocyte': 'Enterocytes',
                'Epithelial cell': 'Epithelial cell',
                'T cell': 'T cell',
                'Stromal cell': 'Stromal',
                'Macrophage': 'Macrophage',
                'B cell': 'B cell',
                'Smooth muscle cell': 'Smooth Muscle',
                'Neutrophil': 'Neutrophil',
                'Endothelial cell (APC)': 'Endothelial',
                'Dendritic cell': 'Dendritic cell',
                'Mast cell': 'Mast cell',
                'Endothelial cell': 'Endothelial',
                'Fetal Neuron': 'Fetal Neuron',
                'Fetal epithelial progenitor': 'Enterocyte Progenitors',
                'Fibroblast': 'Fibroblast',
                'Endothelial cell (endothelial to mesenchymal transition)': 'Endothelial',
                'Fetal stromal cell': 'Stromal',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Monocyte': 'Monocyte',
                'Erythroid cell': 'Erythroid cell',
                'Fetal endocrine cell': 'Enteroendocrine cells',
                'Primordial germ cell': 'Primordial germ cell',
                'Fetal enterocyte': 'Fetal enterocyte',
                'M2 Macrophage': 'Macrophage',
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="AdultTransverseColon_2")
