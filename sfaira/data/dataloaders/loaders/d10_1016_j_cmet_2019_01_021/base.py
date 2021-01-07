import anndata
import numpy as np
import os
import pandas
from typing import Union
from sfaira.data import DatasetBase


class Dataset_d10_1016_j_cmet_2019_01_021(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.download = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117770"

        self.author = "Bhushan"
        self.doi = "10.1016/j.cmet.2019.01.021"
        self.healthy = False
        self.normalization = 'raw'
        self.organ = "pancreas"
        self.organism = "mouse"
        self.protocol = "10x"
        self.state_exact = "diabetic"
        self.sub_tissue = "pancreas"
        self.year = 2019

        self.var_symbol_col = "index"

        self.class_maps = {
            "0": {
                'acinar': 'pancreatic acinar cell',
                'ductal': 'pancreatic ductal cell',
                'leukocyte': 'leukocyte',
                'T cell(Pancreas)': 't cell',
                'B cell(Pancreas)': 'b cell',
                'beta': "pancreatic B cell",
                'alpha': "pancreatic A cell",
                'delta': "pancreatic D cell",
                'pp': "pancreatic PP cell",
                'smooth_muscle': "smooth muscle cell",
                'stellate cell': "pancreatic stellate cell",
                'fibroblast': "stromal cell",
                'endothelial': "endothelial cell"
            },
        }

    def _load_generalized(self, fn, fn_meta):
        celltypes = pandas.read_csv(fn_meta, index_col=0)

        self.adata = anndata.read_mtx(fn + '_matrix.mtx.gz').transpose()
        self.adata.var_names = np.genfromtxt(fn + '_genes.tsv.gz', dtype=str)[:, 1]
        self.adata.obs_names = np.genfromtxt(fn + '_barcodes.tsv.gz', dtype=str)
        self.adata.var_names_make_unique()
        self.adata = self.adata[celltypes.index]
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = celltypes
