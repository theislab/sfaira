import anndata as ad
import os
import pandas as pd
import numpy as np
import scipy


def load(data_dir, sample_fn, **kwargs):
    fn_expr_dict = {
        "11a.6mon": "expression_11a.6mon.txt.gz",
        "GM.6mon": "expression_GM.6mon.txt.gz",
        "HUES66.3mon": "expression_HUES66.3mon.txt",
        "PGP1.3mon.batch2": "expression_PGP1.3mon.batch2.txt.gz",
        "PGP1.3mon": "expression_PGP1.3mon.txt.gz",
        "PGP1.6mon.batch3": "expression_PGP1.6mon.Batch2.txt.gz",
        "PGP1.6mon": "expression_PGP1.6mon.txt.gz",
    }
    fn_tsne_dict = {
        "11a.6mon": "tsne_11a.6mon.1.txt",
        "GM.6mon": "tsne_GM.6mon.1.txt",
        "HUES66.3mon": "tsne_HUES66.3mon.txt",
        "PGP1.3mon.batch2": "tsne_PGP1.3mon.batch2.txt",
        "PGP1.3mon": "tsne_PGP1.3mon.txt",
        "PGP1.6mon.batch3": "tsne_PGP1.6mon.Batch2.1.txt",
        "PGP1.6mon": "tsne_PGP1.6mon.1.txt",
    }
    age_dict = {
        "11a.6mon": "166",
        "GM.6mon": "190",
        "HUES66.3mon": "109",
        "PGP1.3mon.batch2": "113",
        "PGP1.3mon": "101",
        "PGP1.6mon.batch3": "166",
        "PGP1.6mon": "166",
    }

    fn = os.path.join(data_dir, "SCP282")
    fn_meta = os.path.join(fn, "metadata", "meta_combined.txt")
    fn_tsne = os.path.join(fn, "cluster", fn_tsne_dict[sample_fn])
    fn_expression = os.path.join(fn, "expression", fn_expr_dict[sample_fn])

    expr_df = pd.read_csv(fn_expression, sep="\t", index_col=0)
    expr_df = expr_df.apply(np.expm1)
    expr_df = expr_df.T

    metadata_df = pd.read_csv(fn_meta, sep="\t", index_col=0, header=0, low_memory=False).tail(-1)
    metadata_df.index.name = None
    metadata_df = metadata_df.loc[expr_df.index.tolist()].copy()
    metadata_df["organoid_age_days"] = age_dict[sample_fn]

    tsne_df = pd.read_csv(fn_tsne, sep="\t", skiprows=2, index_col=0, header=None)

    # reverse normalisation
    x = scipy.sparse.csr_matrix(np.multiply(expr_df.values, metadata_df["nUMI"].values.reshape(-1, 1).astype(int)), dtype=np.float32)
    x /= 10**6
    x = np.round(x)

    adata = ad.AnnData(
        X=x,
        obs=metadata_df,
        obsm={"X_tsne": tsne_df.loc[metadata_df.index].values},
        var=pd.DataFrame(index=expr_df.columns.values))

    return adata
