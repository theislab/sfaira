import os
import scanpy as sc

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_tar = os.path.join(data_dir, "GSE173351_RAW.tar")
    fn_tar = buffered_decompress(fn_tar)
    fn_tar_inner = os.path.join(fn_tar, sample_fn)
    fn_tar_inner = buffered_decompress(fn_tar_inner)
    # This tar opens up into a secondary folder with a related name:
    fn_unpacked = "_".join(sample_fn.split("_")[1:]).replace(".tar.gz", "")
    adata = sc.read_10x_mtx(os.path.join(fn_tar_inner, fn_unpacked))
    return adata
