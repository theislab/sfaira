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
        DatasetHcl.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.id = "human_trachea_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'trachea'
        self.sub_tissue = 'AdultTrachea'
        self.dev_stage = 'Adult'
        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="AdultTrachea_2")
