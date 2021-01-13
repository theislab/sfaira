import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class.

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
        self.id = "human_rectum_2019_10x_wang_001_10.1084/jem.20191130"

        self.download = "https://covid19.cog.sanger.ac.uk/wang20_rectum.processed.h5ad"

        self.author = "Chen"
        self.doi = "10.1084/jem.20191130"
        self.healthy = True
        self.normalization = 'raw'
        self.organ = "rectum"
        self.organism = "human"
        self.protocol = '10x'
        self.state_exact = 'healthy'
        self.sub_tissue = "rectum"
        self.year = 2019

        self.var_symbol_col = 'index'

        self.obs_key_cellontology_original = 'CellType'

        self.class_maps = {
            "0": {
                'Progenitor': 'Enterocyte progenitor',
                'Goblet': 'Goblet',
                'Enterocyte': 'Enterocyte',
                'Paneth-like': 'Paneth-like',
                'Stem Cell': 'Stem Cell',
                'TA': 'TA',
                'Enteriendocrine': 'Enteroendocrine',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "rectum", "wang20_rectum.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                   .multiply(1 / 10000)
