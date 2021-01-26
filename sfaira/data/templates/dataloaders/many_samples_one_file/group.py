from typing import Union

from sfaira.data import DatasetGroupLoadingOneFile

from .your_dataset_file import Dataset  # ToDo Add your file name here.


class Group(DatasetGroupLoadingOneFile):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None
    ):
        sample_ids = [  # ToDo Add correct sample IDs here.
            "your_sample_id_1",
            "your_sample_id_2"
        ]
        super().__init__(
            cls=Dataset,
            sample_ids=sample_ids,
            path=path, meta_path=meta_path, cache_path=cache_path
        )
