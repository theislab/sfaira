from typing import Union

from sfaira.data import DatasetGroupLoadingManyFiles

from .your_dataset_file import Dataset  # ToDo Add your file name here.


class Group(DatasetGroupLoadingManyFiles):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None
    ):
        sample_fns = [  # ToDo Add correct sample file names here.
            "your_sample_fn_1",
            "your_sample_fn_2"
        ]
        super().__init__(
            cls=Dataset,
            sample_fns=sample_fns,
            path=path, meta_path=meta_path, cache_path=cache_path
        )
