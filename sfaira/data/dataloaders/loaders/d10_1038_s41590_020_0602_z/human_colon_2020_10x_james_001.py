import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class. This dataloader only provides the subset of the published sata which has been made available through the
    covid-19 Cell Atlas.

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
        self.id = "human_colon_2019_10x_james_001_10.1038/s41590-020-0602-z"

        self.download = "https://covid19.cog.sanger.ac.uk/james20.processed.h5ad"
        self.download_meta = None

        self.author = "Teichmann"
        self.doi = "10.1038/s41590-020-0602-z"
        self.healthy = True
        self.normalization = 'raw'
        self.organ = "colon"
        self.organism = "human"
        self.protocol = '10x'
        self.state_exact = 'healthy'
        self.year = 2020

        self.var_symbol_col = 'index'
        self.var_ensembl_col = 'gene_ids'

        self.obs_key_cellontology_original = 'cell_type'

        self.class_maps = {
            "0": {
                'Activated CD4 T': 'Activated CD4 T',
                'B cell IgA Plasma': 'B cell IgA Plasma',
                'B cell IgG Plasma': 'B cell IgG Plasma',
                'B cell cycling': 'B cell cycling',
                'B cell memory': 'B cell memory',
                'CD8 T': 'CD8 T',
                'Follicular B cell': 'Follicular',
                'ILC': 'ILC',
                'LYVE1 Macrophage': 'LYVE1 Macrophage',
                'Lymphoid DC': 'Lymphoid DC',
                'Macrophage': 'Macrophage',
                'Mast': 'Mast cell',
                'Monocyte': 'Monocyte',
                'NK': 'NK',
                'Tcm': 'Tcm',
                'Tfh': 'Tfh',
                'Th1': 'Th1',
                'Th17': 'Th17',
                'Treg': 'Treg',
                'cDC1': 'DC1',
                'cDC2': 'DC2',
                'cycling DCs': 'cycling DCs',
                'cycling gd T': 'cycling gd T',
                'gd T': 'gd T',
                'pDC': 'pDC',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "colon", "james20.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                   .multiply(1 / 10000)
