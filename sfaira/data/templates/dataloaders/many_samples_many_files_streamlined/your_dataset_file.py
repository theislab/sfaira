from typing import Union
import os

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [  # ToDo Add correct sample file names here.
    "your_sample_fn_1",
    "your_sample_fn_2"
]


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def __init__(
            self,
            sample_fn: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(
            sample_fn=sample_fn,
            data_path=data_path,
            meta_path=meta_path,
            cache_path=cache_path,
            **kwargs
        )

        self.id = f"sth_{str(SAMPLE_FNS.index(sample_fn)+1).zfill(3)}_doi"  # ToDo: Index the Dataset ID by the file.
        # ToDo Add you meta data here.

    def _load(self):
        # fn = os.path.join(self.data_dir, self.sample_fn)  # ToDo: add the name of the raw file
        # ToDo: load file fn into self.adata, self.sample_fn represents the current filename.
        pass
