import anndata

from sfaira.unit_tests.mock_data.consts import ASSEMBLY_MOUSE
from sfaira.unit_tests.mock_data.utils import create_adata


def load(data_dir, sample_fn, **kwargs) -> anndata.AnnData:
    ncells = 100
    ngenes = 70
    adata = create_adata(celltypes=["acinar cell", "alpha cell", "beta cell"], ncells=ncells, ngenes=ngenes,
                         assembly=ASSEMBLY_MOUSE)
    return adata
