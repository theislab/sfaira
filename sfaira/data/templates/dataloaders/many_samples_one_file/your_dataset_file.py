import anndata
from typing import Union

from sfaira.data import DatasetBaseGroupLoadingManyFiles


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
        self.obs_key_sample = 'Sample'  # ToDo: Make sure to include this attribute which indicates the column in
        # self.adata in which you saved the sample IDs based on which the full adata object is subsetted.

    def _load_full_group_object(self, fn=None) -> anndata.AnnData:
        pass  # ToDo: load full data object and return (no subsetting!)
