import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.io
import os
import gzip


def load(data_dir, sample_fn, **kwargs):
    age_dict = {
        "1M": "30",
        "3M": "90",
        "6M": "180",
        "10M": "300",
    }
    fn = os.path.join(data_dir, "GSE130238_RAW.tar")
    with tarfile.open(fn) as tar:
        with gzip.open(tar.extractfile(f"{sample_fn}_cortical_organoids_matrix.mtx.gz"), "rb") as mm:
            x = scipy.io.mmread(mm).T.tocsr().astype(np.float32)
        obs = pd.read_csv(tar.extractfile(f"{sample_fn}_cortical_organoids_barcodes.tsv.gz"),
                          delimiter='\t', compression="gzip", header=None, index_col=0, names=[None])
        obs['organoid_age_days'] = age_dict[sample_fn.split("_")[1]]
        var = pd.read_csv(tar.extractfile(f"{sample_fn}_cortical_organoids_genes.tsv.gz"),
                          delimiter='\t', compression='gzip', index_col=0, header=None, names=[None, "gene_id"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    return adata
