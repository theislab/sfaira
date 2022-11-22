import os
import scanpy as sc

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    data_dir = buffered_decompress(os.path.join(data_dir, "GSE156989_RAW.tar"))
    fn = os.path.join(data_dir, sample_fn + "_labeled.h5")
    adata = sc.read_10x_h5(fn)
    return adata
