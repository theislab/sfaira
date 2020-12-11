import anndata
import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_TMS


class DatasetTms(DatasetBase):
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
        self._ADATA_IDS_TMS = ADATA_IDS_TMS()
        self.fn = fn

        self.obs_key_cellontology_class = self._ADATA_IDS_TMS.cell_ontology_class
        self.obs_key_cellontology_id = self._ADATA_IDS_TMS.cell_ontology_id
        self.obs_key_cellontology_original = self._ADATA_IDS_TMS.cell_types_original
        self.obs_key_dev_stage = self._ADATA_IDS_TMS.dev_stage
        self.obs_key_state_exact = self._ADATA_IDS_TMS.state_exact

        self.author = "Quake"
        self.year = "2019"
        self.doi = "10.1101/661728"
        self.protocol = None  # TODO load from data / name
        self.normalization = 'norm'
        self.healthy = True
        self.state_exact = "healthy"
        self.species = "mouse"

        self.var_ensembl_col = self._ADATA_IDS_TMS.gene_id_ensembl
        self.var_symbol_col = self._ADATA_IDS_TMS.gene_id_names

    def _load_tms(self, fn):
        self.adata = anndata.read_h5ad(fn)
        if self.source == "aws":
            self.adata.X = self.adata.raw.X
            self.adata.var = self.adata.raw.var
            del self.adata.raw
            self.adata.obsm = {}
            self.adata.varm = {}
            self.adata.uns = {}
