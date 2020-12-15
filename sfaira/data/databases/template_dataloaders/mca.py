import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase


class DatasetMca(DatasetBase):
    """
    This is a dataloader template for mca data.
    """

    def __init__(
            self,
            path: Union[str, None],
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)

        self.obs_key_cellontology_class = "Annotation"
        self.obs_key_cellontology_original = "Annotation"

        self.author = "Guo"
        self.year = "2018"
        self.doi = "10.1016/j.cell.2018.02.001"
        self.protocol = "microwell-seq"
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = "healthy"
        self.species = "mouse"

        self.var_ensembl_col = "ensembl"
        self.var_symbol_col = "names"

    def _load_mca(self, fn, fn_meta):
        celltypes = pandas.read_csv(fn_meta, index_col=1)
        celltypes = celltypes.drop(['Unnamed: 0'], axis=1)

        data = pandas.read_csv(fn, sep=' ', header=0)
        self.adata = anndata.AnnData(data.T)
        self.adata = self.adata[np.array([x in celltypes.index for x in self.adata.obs_names])].copy()
        self.adata.obs = celltypes.loc[self.adata.obs_names, :]
