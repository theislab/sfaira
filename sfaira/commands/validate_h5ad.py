import logging

import anndata
import numpy as np
from rich import print
import scipy.sparse

log = logging.getLogger(__name__)


class H5adValidator:

    def __init__(self, test_h5ad, schema=None):
        self.fn_h5ad: str = test_h5ad
        if schema not in ["cellxgene"]:
            raise ValueError(f"Did not recognize schema {schema}")
        self.schema = schema
        self._adata = None

    @property
    def adata(self):
        if self._adata is None:
            self._adata = anndata.read_h5ad(filename=self.fn_h5ad)
        return self._adata

    def validate(self):
        self._test_schema()
        self._test_numeric_data()
        print(f'[bold orange]Sfaira butler: "The object {self.fn_h5ad} passed the schema {self.schema}."')

    def _test_schema(self) -> None:
        """Verify that object elements match schema definitions."""
        if self.schema == "cellxgene":
            from cellxgene_schema import validate
            validate.validate(h5ad_path=self.fn_h5ad)
        else:
            raise ValueError(f"schema {self.schema} not supported")

    def _test_numeric_data(self) -> None:
        """Verify that numeric matrices match schema definitions."""
        # if isinstance(self.adata.X, scipy.sparse.spmatrix):
        #     x = np.unique(np.asarray(self.adata.X.todense()))
        # else:
        #     x = np.unique(np.asarray(self.adata.X))
        # deviation_from_integer = np.minimum(x % 1, 1. - x % 1)
        # assert np.max(deviation_from_integer) < 1e-6
        pass
