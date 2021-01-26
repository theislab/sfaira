from typing import Union

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [  # ToDo Add correct sample file names here.
    "your_sample_fn_1",
    "your_sample_fn_2"
]


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def __init__(
            self,
            sample_fn: str,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(
            sample_fn=sample_fn,
            path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)

        # ToDo Add you meta data here.

    def _load_any_object(self, fn):
        pass  # ToDo: load file fn into self.adata.
