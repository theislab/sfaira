import os
import scanpy as sc

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_tar = os.path.join(data_dir, "GSE160641_RAW.tar")
    fn_tar = buffered_decompress(fn_tar)
    adata = sc.read_10x_mtx(fn_tar, prefix=sample_fn + "_")
    return adata
