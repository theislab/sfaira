import tarfile
import os


import pandas as pd
import numpy as np
import anndata as ad
from scipy.io import mmread


def load(data_dir, **kwargs):
    dataset_tar_path = os.path.join(data_dir, "GSE123722_RAW.tar")
    with tarfile.open(dataset_tar_path, "r") as tar:
        tar.extractall(data_dir)

    matrix_path = os.path.join(data_dir, "GSM4589792_hSpS_d53_matrix.mtx.gz")
    genes_path = os.path.join(data_dir, "GSM4589792_hSpS_d53_features.tsv.gz")
    barcodes_path = os.path.join(data_dir, "GSM4589792_hSpS_d53_barcodes.tsv.gz")

    genes_df = pd.read_csv(genes_path, compression='gzip', sep="\t", header=None)
    genes_df = genes_df.drop(columns=2)
    genes_df.columns = ["Ensembl", "Gene_Symbol"]

    cells_df = pd.read_csv(barcodes_path, compression='gzip', sep="\t", header=None)
    cells_df.columns = ["Barcode"]

    expression_mt = mmread(matrix_path)

    adata = ad.AnnData(X=expression_mt.T, obs=cells_df, var=genes_df, dtype=np.int64)

    return adata
