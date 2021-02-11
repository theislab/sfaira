from typing import Union
import os
from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)

        # ToDo: Add you meta data here.

    def _load(self):
        fn = os.path.join(self.doi_path, )  # ToDo: add the name of the raw file
        # ToDo: add code that loads to raw file into an AnnData object
