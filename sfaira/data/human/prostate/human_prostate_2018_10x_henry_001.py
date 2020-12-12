import anndata
import os
from typing import Union
from .external import DatasetBase
import numpy as np
import scipy.sparse


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
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_prostate_2018_10x_henry_001_10.1016/j.celrep.2018.11.086"
        self.download = "https://covid19.cog.sanger.ac.uk/henry18_0.processed.h5ad"
        self.download_meta = None
        self.organ = "prostate"
        self.sub_tissue = "prostate"
        self.author = "Strand"
        self.year = 2018
        self.doi = "10.1016/j.celrep.2018.11.086"
        self.protocol = '10x'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'CellType'

        self.class_maps = {
            "0": {
                'Basal': 'Basal cell',
                'Hillock': 'Hillock',
                'Luminal': 'Luminal',
                'Endothelia': 'Endothelial cell',
                'Club': 'Club',
                'Fibroblast': 'Fibroblast',
                'Smooth muscle': 'Smooth muscle cell',
                'Leukocytes': 'Leukocytes',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "prostate", "henry18_0.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)
