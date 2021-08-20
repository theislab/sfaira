import anndata

from sfaira.unit_tests.data_for_tests.loaders.consts import ASSEMBLY_HUMAN, CELLTYPES
from sfaira.unit_tests.data_for_tests.loaders.utils import _create_adata


def load(data_dir, sample_fn, **kwargs) -> anndata.AnnData:
    ncells = 100
    ngenes = 60
    adata = _create_adata(celltypes=CELLTYPES[:2], ncells=ncells, ngenes=ngenes,
                          assembly=ASSEMBLY_HUMAN)
    return adata
