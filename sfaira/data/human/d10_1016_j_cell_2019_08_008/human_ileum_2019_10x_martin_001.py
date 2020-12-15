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
        super().__init__(path=path, meta_path=meta_path, **kwargs)
        self.organism = "human"
        self.id = "human_ileum_2019_10x_martin_001_10.1016/j.cell.2019.08.008"
        self.download = "https://covid19.cog.sanger.ac.uk/martin19.processed.h5ad"
        self.download_meta = None
        self.organ = "ileum"
        self.sub_tissue = "ileum"
        self.author = "Kenigsberg"
        self.year = 2019
        self.doi = "v"
        self.protocol = '10x'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.var_ensembl_col = 'gene_ids'
        self.obs_key_cellontology_original = 'CellType'

        self.class_maps = {
            "0": {
                'T cells': 'T cells',
                'Plasma cells': 'Plasma Cells',
                'B cells': 'B cells',
                'MNP': 'MNP',
                'ILC': 'ILC',
                'Enterocytes': 'Enterocytes',
                'Fibs': 'Fibroblasts',
                'CD36+ endothelium': 'CD36+ endothelium',
                'Progenitors': 'Progenitors',
                'Goblets': 'Goblet cells',
                'Glial cells': 'Glial cells',
                'Cycling': 'Cycling',
                'ACKR1+ endothelium': 'ACKR1+ endothelium',
                'Pericytes': 'Pericytes',
                'Lymphatics': 'Lymphatics',
                'Mast cells': 'Mast cells',
                'SM': 'Smooth muscle cell',
                'TA': 'TA',
                'Paneth cells': 'Paneth cells',
                'Enteroendocrines': 'Enteroendocrine cells',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "ileum", "martin19.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                   .multiply(1/10000)
        self.adata = self.adata[self.adata.obs['CellType'] != 'Doublets'].copy()
