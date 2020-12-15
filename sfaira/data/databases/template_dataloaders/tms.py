import anndata
import os
from typing import Union
from .external import DatasetBase


class DatasetTms(DatasetBase):
    """
    This is a dataloader template for tabula muris data.
    """

    def __init__(
            self,
            path: Union[str, None],
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)

        self.obs_key_cellontology_class = "cell_ontology_class"
        self.obs_key_cellontology_id = "cell_ontology_id"
        self.obs_key_cellontology_original = "cell_types_original"
        self.obs_key_age = "age"
        self.obs_key_dev_stage = "development_stage"

        self.author = "Quake"
        self.year = "2019"
        self.doi = "10.1101/661728"
        self.normalization = 'norm'
        self.healthy = True
        self.state_exact = "healthy"
        self.organism = "mouse"

        self.var_ensembl_col = None
        self.var_symbol_col = "index"

    def _load_tms(self, fn):
        self.adata = anndata.read_h5ad(fn)
        if self.source == "aws":
            self.adata.X = self.adata.raw.X
            self.adata.var = self.adata.raw.var
            del self.adata.raw
            self.adata.obsm = {}
            self.adata.varm = {}
            self.adata.uns = {}

    def _get_protocol_tms(self, x) -> str:
        return "smartseq2" if "smartseq2" in x else "10x"
