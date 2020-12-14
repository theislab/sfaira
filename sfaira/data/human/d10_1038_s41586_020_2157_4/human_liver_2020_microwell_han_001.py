import os
from typing import Union
from .external import DatasetHcl


class Dataset(DatasetHcl):
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
        DatasetHcl.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.id = "human_liver_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'Liver'
        self.sub_tissue = 'AdultLiver'
        self.dev_stage = 'Adult'
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
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "liver", "hcl_AdultLiver_1.h5ad")
            self._load_hcl(fn=fn)
