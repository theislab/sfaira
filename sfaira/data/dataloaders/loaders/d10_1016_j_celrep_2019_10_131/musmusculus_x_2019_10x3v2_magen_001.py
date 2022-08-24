import anndata
import os
import pandas as pd
import scipy.io

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_features = os.path.join(data_dir, "GSE124691_Genes.tsv.gz")
    fn_tar = os.path.join(data_dir, "GSE124691_RAW.tar")
    fn_tar = buffered_decompress(fn_tar)
    # fn_bc = os.path.join(fn_tar, sample_fn + "_barcodes.tsv.gz")
    fn_x = os.path.join(fn_tar, sample_fn + "_matrix.mtx.gz")

    x = scipy.io.mmread(fn_x).T.tocsr()
    # obs = pd.read_csv(fn_bc, header=None, sep="\t", index_col=0, compression="gzip")
    var = pd.read_csv(fn_features, header=None, sep="\t", index_col=0, names=['symbols'], compression="gzip")
    adata = anndata.AnnData(X=x, var=var)
    # Causes error? and not used.
    # adata.obs_names = obs.index
    return adata
