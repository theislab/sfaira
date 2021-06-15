import anndata

from sfaira.unit_tests.mock_data.consts import ASSEMBLY_HUMAN
from sfaira.unit_tests.mock_data.utils import create_adata


def load(data_dir, sample_fn, **kwargs) -> anndata.AnnData:
    ncells = 100
    ngenes = 60
    adata = create_adata(celltypes=["adventitial cell", "endothelial cell"], ncells=ncells, ngenes=ngenes,
                         assembly=ASSEMBLY_HUMAN)
    return adata
