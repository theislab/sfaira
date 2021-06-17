import logging

import anndata
import numpy as np
import scipy.sparse

log = logging.getLogger(__name__)


class H5adValidator:

    def __init__(self, test_h5ad, schema=None):
        self.fn: str = test_h5ad
        if schema not in ["cellxgene"]:
            raise ValueError(f"Did not recognize schema {schema}")
        self.schema = schema
        self._adata = None

    @property
    def adata(self):
        if self.adata is None:
            self._adata = anndata.read_h5ad(filename=self.fn)
        return self._adata

    def test_schema(self) -> None:
        """Verify that object elements match schema definitions."""
        if self.schema == "cellxgene":
            from cellxgene.backend.server.converters.schema.validate import validate_adata
            validate_adata(adata=self.adata, shallow=False)
        else:
            assert False

    def test_numeric_data(self) -> None:
        """Verify that numeric matrices match schema definitions."""
        if isinstance(self.adata.X, scipy.sparse.spmatrix):
            x = np.unique(np.asarray(self.adata.X.todense()))
        else:
            x = np.unique(np.asarray(self.adata.X))
        deviation_from_integer = np.minimum(x % 1, 1. - x % 1)
        print(f"Maximum deviation from integer found: {np.max(deviation_from_integer)}")
        assert np.max(deviation_from_integer) < 1e-3
