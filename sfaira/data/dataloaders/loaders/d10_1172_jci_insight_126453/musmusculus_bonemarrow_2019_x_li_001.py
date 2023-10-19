import anndata
import os
import pandas as pd
import scanpy as sc
import scipy.io

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_tar = os.path.join(data_dir, "GSE117176_RAW.tar")
    fn_tar = buffered_decompress(fn_tar)
    # Note: right now, scanpy.read_10x_mtx does not discover this is an old data format with a
    # genes.tsv instead of features.tsv because this file has a gz postfix
    # Therefore need to decompress first:
    _ = buffered_decompress(os.path.join(fn_tar, sample_fn + "_barcodes.tsv.gz"))
    fn_genes = buffered_decompress(os.path.join(fn_tar, sample_fn + "_genes.tsv.gz"))
    # One matrix file has a spelling error in the name:
    if sample_fn == "GSM3272969_M1_BMDM":
        sample_fn_mtx = "GSM3272969_M1_BDMD"
        fn_mtx = buffered_decompress(os.path.join(fn_tar, sample_fn_mtx + "_matrix.mtx.gz"))
        x = scipy.io.mmread(fn_mtx).T.tocsr()
        var = pd.read_csv(fn_genes, header=None, sep="\t", names=['gene_ids', 'symbols'])
        var.index = var['symbols'].values
        del var['symbols']
        adata = anndata.AnnData(X=x, var=var)
    else:
        _ = buffered_decompress(os.path.join(fn_tar, sample_fn + "_matrix.mtx.gz"))
        adata = sc.read_10x_mtx(fn_tar, prefix=sample_fn + "_")
    return adata
