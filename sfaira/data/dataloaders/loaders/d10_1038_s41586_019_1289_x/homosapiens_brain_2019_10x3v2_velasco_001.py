import anndata as ad
import os
import pandas as pd
import numpy as np
import scipy


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    fn_meta = os.path.join(data_dir, "meta_combined.txt")

    expr_df = pd.read_csv(fn, sep="\t", compression='gzip', index_col=0)
    expr_df = expr_df.apply(np.expm1)
    expr_df = expr_df.T

    metadata_df = pd.read_csv(fn_meta, sep="\t", index_col=0, header=0)
    metadata_df = metadata_df.tail(-1)

    # multiply by the UMIs, divide by 10^6
    arr = scipy.sparse.csr_matrix(
        np.multiply(expr_df.values, metadata_df.nUMI.loc[
            expr_df.index.tolist()].values.reshape(-1, 1).astype(int)))
    arr /= 10**6

    adata = ad.AnnData(
        X=arr, obs=metadata_df.loc[
            expr_df.index.tolist()].copy(),
        var=pd.DataFrame(index=expr_df.columns.values))

    return adata
