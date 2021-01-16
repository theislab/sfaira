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
        self.id = "human_thymus_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'thymus'
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
        self._load_generalized(fn=fn, sample_id="FetalThymus_2")
