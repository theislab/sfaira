import os
import pandas as pd
import scanpy as sc

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    # Note: right now, scanpy.read_10x_mtx does not discover this is an old data format with a
    # genes.tsv instead of features.tsv because this file has a gz postfix
    # Therefore need to decompress first:
    _ = buffered_decompress(os.path.join(data_dir, sample_fn + "_barcodes.tsv.gz"))
    _ = buffered_decompress(os.path.join(data_dir, sample_fn + "_genes.tsv.gz"))
    _ = buffered_decompress(os.path.join(data_dir, sample_fn + "_matrix.mtx.gz"))
    if sample_fn == "GSE114687_HuMEC":
        fn_meta = os.path.join(data_dir, sample_fn + "_pseudospace_metadata.tsv.gz")
    else:
        fn_meta = os.path.join(data_dir, sample_fn + "_metadata.tsv.gz")
    adata = sc.readwrite.read_10x_mtx(data_dir, prefix=sample_fn + "_")
    # Meta data files all look different, give all the same column space:
    if sample_fn == "GSE114687_CROPSeq_pseudospace":
        tab = pd.read_csv(fn_meta, header=None, index_col=None, sep="\t")
        tab["location"] = tab.iloc[:, 15].values
        adata.obs = tab
    elif sample_fn == "GSE114687_HuMEC":
        tab = pd.read_csv(fn_meta, header=0, index_col=0, sep="\t")
        tab["location"] = tab["spatial_id"].values
        adata.obs = tab
    elif sample_fn == "GSE114687_pseudospace":
        pass
    else:
        assert False
    return adata
