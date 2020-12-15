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
        self.id = "human_thymus_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'thymus'
        self.sub_tissue = 'FetalThymus'
        self.dev_stage = 'Fetus'
        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'B cell': 'B cell',
                'CB CD34+': 'CB CD34+',
                'Dendritic cell': 'Dendritic cell',
                'Erythroid cell': 'Ery',
                'Erythroid progenitor cell (RP high)': 'Ery',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Macrophage': 'Mac',
                'Monocyte': 'Mono',
                'Neutrophil': 'Neutrophil',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Proliferating T cell': 'Proliferating T cell',
                'T cell': 'T cell',
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="FetalThymus_2")
