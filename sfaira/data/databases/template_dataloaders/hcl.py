import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_HCL


class DatasetHcl(DatasetBase):
    """
    This is a dataloader template for tabula muris data.
    """

    def __init__(
            self,
            path: Union[str, None],
            fn: str,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self._ADATA_IDS_HCL = ADATA_IDS_HCL()
        self.fn = fn

        self.obs_key_cellontology_class = self._ADATA_IDS_HCL.cell_ontology_class
        self.obs_key_cellontology_original = self._ADATA_IDS_HCL.cell_types_original

        self.author = 'Guo'
        self.doi = '10.1038/s41586-020-2157-4'
        self.species = "human"
        self.download = 'https://figshare.com/articles/HCL_DGE_Data/7235471'
        self.download_meta = None
        self.healthy = True
        self.normalization = 'raw'
        self.protocol = 'microwell-seq'
        self.state_exact = 'healthy'
        self.year = 2020

        self.var_ensembl_col = self._ADATA_IDS_HCL.gene_id_ensembl
        self.var_symbol_col = self._ADATA_IDS_HCL.gene_id_names

    def _load_hcl(self, fn):
        pass  # TODO
