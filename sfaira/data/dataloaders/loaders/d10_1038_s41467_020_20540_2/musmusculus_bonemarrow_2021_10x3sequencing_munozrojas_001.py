import os
import scanpy as sc

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_tar = os.path.join(data_dir, "GSE161125_RAW.tar")
    fn_tar = buffered_decompress(fn_tar)
    # Note: right now, scanpy.read_10x_mtx does not discover this is an old data format with a
    # genes.tsv instead of features.tsv because this file has a gz postfix
    # Therefore need to decompress first:
    _ = buffered_decompress(os.path.join(fn_tar, sample_fn + "_barcodes.tsv.gz"))
    _ = buffered_decompress(os.path.join(fn_tar, sample_fn + "_genes.tsv.gz"))
    _ = buffered_decompress(os.path.join(fn_tar, sample_fn + "_matrix.mtx.gz"))
    adata = sc.read_10x_mtx(fn_tar, prefix=sample_fn + "_")
    return adata
