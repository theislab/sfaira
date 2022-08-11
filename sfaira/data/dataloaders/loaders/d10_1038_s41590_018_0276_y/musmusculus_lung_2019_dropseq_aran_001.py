import anndata
import os

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    dir_tar = os.path.join(data_dir, "GSE111664_RAW.tar")
    dir_tar = buffered_decompress(dir_tar)
    fn = os.path.join(dir_tar, sample_fn)
    # replace the simple data loading code below with the code required to load your data file(s)
    adata = anndata.read_text(fn).T

    return adata
