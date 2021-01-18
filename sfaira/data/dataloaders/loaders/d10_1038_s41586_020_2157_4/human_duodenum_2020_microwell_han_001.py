from typing import Union
from .base import Dataset_d10_1038_s41586_020_2157_4


class Dataset(Dataset_d10_1038_s41586_020_2157_4):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_duodenum_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = "duodenum"
        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        self._load_generalized(fn=fn, sample_id="AdultDuodenum_1")
