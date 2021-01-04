import os
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
        self.id = "human_omentum_2020_microwell_han_002_10.1038/s41586-020-2157-4"
        self.organ = 'omentum'
        self.sub_tissue = 'AdultOmentum'
        self.dev_stage = 'Adult'
        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, sample_id="AdultOmentum_3")
