import os
import scanpy as sc

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_tar = os.path.join(data_dir, "GSE143317_RAW.tar")
    fn_tar = buffered_decompress(fn_tar)
    fn_tar_inner = os.path.join(fn_tar, f'{sample_fn}_filtered_feature_bc_matrix.tar.gz')
    fn_tar_inner = buffered_decompress(fn_tar_inner)
    adata = sc.read_10x_mtx(fn_tar_inner, var_names='gene_ids')
    return adata
