import anndata
import os

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    dir_tar = os.path.join(data_dir, "GSE176031_RAW.tar")
    dir_tar = buffered_decompress(dir_tar)
    fn = os.path.join(dir_tar, sample_fn)
    adata = anndata.read_text(fn).T

    return adata
