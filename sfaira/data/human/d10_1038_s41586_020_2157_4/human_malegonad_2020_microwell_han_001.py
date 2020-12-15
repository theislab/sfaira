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
        self.id = "human_malegonad_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'malegonad'
        self.sub_tissue = 'FetalMaleGonad'
        self.dev_stage = 'Fetus'
        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'B cell': 'B cell',
                'CB CD34+': 'CB CD34+',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cells',
                'Erythroid cell': 'Erythroid cell',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Fasciculata cell': 'Fasciculata cell',
                'Fetal acinar cell': 'Fetal acinar cell',
                'Fetal chondrocyte': 'Fetal chondrocyte',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Fetal fibroblast': 'Fetal fibroblast',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Fetal neuron': 'Fetal neuron',
                'Fetal skeletal muscle cell': 'Fetal skeletal muscle cell',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Immature sertoli cell (Pre-Sertoli cell)': 'Sertoli cells',
                'Loop of Henle': 'Loop of Henle',
                'Macrophage': 'Macrophages',
                'Monocyte': 'Monocyte',
                'Neutrophil': 'Neutrophil',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Primordial germ cell': 'Primordial germ cell',
                'Proximal tubule progenitor': 'Proximal tubule progenitor',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
                'Ureteric bud cell': 'Ureteric bud cell',
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="FetalMaleGonad_1")
