import os
import scanpy as sc

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, **kwargs):
    # Note: right now, scanpy.read_10x_mtx does not discover this is an old data format with a
    # genes.tsv instead of features.tsv because this file has a gz postfix
    # Therefore need to decompress first:
    _ = buffered_decompress(os.path.join(data_dir, "barcodes.tsv.gz"))
    _ = buffered_decompress(os.path.join(data_dir, "genes.tsv.gz"))
    _ = buffered_decompress(os.path.join(data_dir, "matrix.mtx.gz"))
    adata = sc.read_10x_mtx(data_dir)
    del adata.var["gene_ids"]
    return adata
