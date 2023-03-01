import anndata as ad
import os
import scipy.io
import pandas as pd
import numpy as np
import tarfile
import gzip


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "rna_matrices.tar.gz")
    with tarfile.open(fn) as tar:
        with gzip.open(tar.extractfile("data_matrices/counts.mtx.gz"), "rb") as mm:
            x = scipy.io.mmread(mm).T.tocsr().astype(np.float32)
        var = pd.read_csv(tar.extractfile("data_matrices/features.tsv.gz"), compression="gzip", sep="\t", header=None, index_col=0, names=[None])
        obs = pd.read_csv(tar.extractfile("data_matrices/barcodes.tsv.gz"), compression="gzip", sep="\t", index_col=0, header=None, names=[None])
        meta = pd.read_csv(tar.extractfile("data_matrices/meta.tsv.gz"), compression="gzip", sep="\t", index_col=0)

    meta = meta.loc[obs.index]
    meta["organoid_age_days"] = [i.split("_")[1] for i in meta["sample"]]
    meta["celltype"] = meta["lineage"] + "-" + meta["stage_manual"]

    adata = ad.AnnData(X=x, obs=meta, var=var)

    return adata
