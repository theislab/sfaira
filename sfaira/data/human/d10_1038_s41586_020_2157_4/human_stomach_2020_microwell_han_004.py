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
        self.id = "human_stomach_2020_microwell_han_004_10.1038/s41586-020-2157-4"
        self.organ = 'stomach'
        self.sub_tissue = 'FetalIntestine'
        self.dev_stage = 'Fetus'
        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="FetalIntestine_3")
