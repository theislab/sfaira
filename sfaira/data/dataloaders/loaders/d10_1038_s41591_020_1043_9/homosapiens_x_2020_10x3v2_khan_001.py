import os
import numpy as np
import pandas as pd
import gzip
import scipy
import tarfile
import anndata as ad
import scipy.io


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE145122_RAW.tar")
    with tarfile.open(fn) as tar:
        obs = pd.read_csv(tar.extractfile(f"{sample_fn}_barcodes.tsv.gz"),
                          compression="gzip", header=None, index_col=0, names=[None])
        var = pd.read_csv(tar.extractfile(f"{sample_fn}_features.tsv.gz"),
                          compression="gzip", sep='\t', header=None, usecols=[0, 1], index_col=1, names=["ensembl_id", None])
        with gzip.open(tar.extractfile(f"{sample_fn}_matrix.mtx.gz"), 'rb') as mm:
            x = scipy.io.mmread(mm).T.tocsr().astype(np.float32)
    obs["organoid_age_days"] = "83"
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata
