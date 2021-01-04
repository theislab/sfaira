from typing import Union
from .base import Dataset_d10_1038_s41586_020_2157_4


class Dataset(Dataset_d10_1038_s41586_020_2157_4):
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
        self.id = "human_liver_2020_microwell_han_005_10.1038/s41586-020-2157-4"
        self.organ = 'Liver'
        self.sub_tissue = 'FetalLiver'
        self.dev_stage = 'Fetus'
        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'B cell': 'B cell',
                'B cell (Plasmocyte)': 'Plasma B cell',
                'CB CD34+': 'CB CD34+',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cell',
                'Endothelial cell (APC)': 'Endothelial cell (APC)',
                'Endothelial cell (endothelial to mesenchymal transition)': 'Endothelial cell (endothelial to mesenchymal transition)',
                'Enterocyte progenitor': 'Enterocyte progenitor',
                'Erythroid cell': 'Late Erythroid',
                'Erythroid progenitor cell (RP high)': 'Early Erythroid',
                'Fetal enterocyte ': 'Enterocyte ',
                'Fetal epithelial progenitor': 'Epithelial progenitor',
                'Fetal fibroblast': 'Fibroblast',
                'Gastric endocrine cell': 'Gastric endocrine cell',
                'Goblet cell': 'Goblet cell',
                'Macrophage': 'Non inflammatory macrophages',
                'Mast cell': 'Mast cell',
                'Monocyte': 'Monocyte',
                'Myeloid cell': 'Myeloid cell',
                'Neutrophil': 'Neutrophil',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Pancreas exocrine cell': 'Pancreas exocrine cell',
                'Primordial germ cell': 'Primordial germ cell',
                'Proliferating T cell': 'Proliferating T cell',
                'Sinusoidal endothelial cell': 'Liver sinusoidal endothelial cells',
                'Smooth muscle cell': 'Smooth muscle cell',
                'T cell': 'T cell'
            },
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, sample_id="Liver_2")
