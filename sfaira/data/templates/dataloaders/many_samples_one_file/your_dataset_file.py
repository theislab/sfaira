import anndata
from typing import Union
import os

from sfaira.data import DatasetBaseGroupLoadingOneFile

SAMPLE_IDS = [  # ToDo Add correct sample IDs here.
    "your_sample_id_1",
    "your_sample_id_2"
]


class Dataset(DatasetBaseGroupLoadingOneFile):

    def __init__(
            self,
            sample_id: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(
            sample_id=sample_id,
            data_path=data_path,
            meta_path=meta_path,
            cache_path=cache_path,
            **kwargs
        )

        # ToDo Add you meta data here.
        self.obs_key_sample = 'Sample'  # ToDo: Make sure to include this attribute which indicates the column in
        # self.adata in which you saved the sample IDs based on which the full adata object is subsetted.

    def _load_full(self) -> anndata.AnnData:
        fn = os.path.join(self.data_dir, )  # ToDo: add the name of the raw file
        # ToDo: load full data into AnnData object (no subsetting!)
