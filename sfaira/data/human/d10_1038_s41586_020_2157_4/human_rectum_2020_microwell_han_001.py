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
        self.id = "human_rectum_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'rectum'
        self.sub_tissue = 'AdultRectum'
        self.dev_stage = 'Adult'
        self.class_maps = {
            "0": {
                'B cell': 'B cell',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell (APC)': 'Endothelial cell (APC)',
                'Enterocyte': 'Enterocyte',
                'Enterocyte progenitor': 'Enterocyte progenitor',
                'Epithelial cell': 'Epithelial cell',
                'Erythroid cell': 'Erythroid cell',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Macrophage': 'Macrophage',
                'Mast cell': 'Mast cell',
                'Monocyte': 'Monocyte',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="AdultRectum_1")
