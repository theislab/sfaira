import anndata
import os
from typing import Union

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.download_url_data = "private,fetal_liver_alladata_.h5ad"
        self.download_url_meta = None

        self.author = "Popescu"
        self.doi = "10.1038/s41586-019-1652-y"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "liver"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"
        self.obs_key_cellontology_original = "cell.labels"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "fetal_liver_alladata_.h5ad")
    adata = anndata.read(fn)

    return adata
