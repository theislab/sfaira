import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.author = "Bhushan"
        self.id = "mouse_pancreas_2019_10x_thompson_004_10.1016/j.cmet.2019.01.021"
        self.doi = "10.1016/j.cmet.2019.01.021"
        self.download = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117770"
        self.healthy = False
        self.normalization = 'raw'
        self.organ = "pancreas"
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

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "pancreas", "GSM3308549_NOD_14w_B")
            fn_meta = os.path.join(self.path, "mouse", "pancreas", "GSM3308549_NOD_14w_B_annotation.csv")

        celltypes = pandas.read_csv(fn_meta, index_col=0)

        self.adata = anndata.read_mtx(fn + '_matrix.mtx.gz').transpose()
        self.adata.var_names = np.genfromtxt(fn + '_genes.tsv.gz', dtype=str)[:, 1]
        self.adata.obs_names = np.genfromtxt(fn + '_barcodes.tsv.gz', dtype=str)
        self.adata.var_names_make_unique()
        self.adata = self.adata[celltypes.index]
        self.set_unkown_class_id(ids=[np.nan, "nan"])
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = celltypes
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = celltypes
