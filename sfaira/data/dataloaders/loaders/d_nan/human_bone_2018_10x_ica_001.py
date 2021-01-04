import anndata
import os
from typing import Union
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader supports reading of the downloaded raw data file if `load_raw=True` is passed to self.load()
    To download the datafile required by this dataloader, use the link provided as the `download_website` attribute of
    this class. For (up to 100-fold faster) repeated data loading, please pass `load_raw=False` when calling the
    self.load() method. For this, you need to preprocess the raw files as below and place the resulting h5ad file in the
    data folder of this organ:

    import anndata
    import numpy as np
    adata = anndata.read_loom('c95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom')
    idx = np.logical_and((adata.obs['derived_organ_parts_label'] == 'bone marrow').values, (adata.obs['emptydrops_is_cell'] == 't').values)
    adata = adata[idx].copy()
    adata.write('ica_bone.h5ad')

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):

        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.organism = "loaders"
        self.id = "human_bone_2018_10x_ica_unknown"
        self.download = "https://data.humancellatlas.org/project-assets/project-matrices/cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom"
        self.download_meta = None
        self.organ = "bone"
        self.sub_tissue = "bone_marrow"
        self.author = 'Regev'
        self.year = 2018
        self.doi = None
        self.protocol = '10x'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.var_ensembl_col = 'Accession'

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "loaders", "bone", "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom")
        self.adata = anndata.read_loom(fn)
        idx = np.logical_and((self.adata.obs['derived_organ_parts_label'] == 'bone marrow').values,
                             (self.adata.obs['emptydrops_is_cell'] == 't').values)
        self.adata = self.adata[idx].copy()
