import tarfile
import os
import pandas as pd
import numpy as np
import anndata as ad
import scipy.io
import scipy.sparse
import gzip


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE123722_RAW.tar")
    with tarfile.open(fn, "r") as tar:
        if sample_fn == "GSM4589792_hSpS_d53":
            obs = pd.read_csv(tar.extractfile(f"{sample_fn}_barcodes.tsv.gz"),
                              compression="gzip", header=None, index_col=0, names=[None])
            var = pd.read_csv(tar.extractfile(f"{sample_fn}_features.tsv.gz"),
                              compression="gzip", sep='\t', header=None, usecols=[0, 1], index_col=1, names=["ensembl_id", None])
            with gzip.open(tar.extractfile(f"{sample_fn}_matrix.mtx.gz"), 'rb') as mm:
                x = scipy.io.mmread(mm).T.tocsr().astype(np.float32)
            adata = ad.AnnData(X=x, obs=obs, var=var)
            adata.var_names_make_unique()
            adata.obs["cellline"] = adata.obs.index.str.split("-").str[1]
            adata.obs["cellline"].replace({"1": "8858-1", "2": "1205-4", "3": "0524-1"}, inplace=True)
        else:
            df = pd.read_csv(tar.extractfile(f"{sample_fn}_MolsPerCell.csv.gz"), compression="gzip", index_col=0)
            adata = ad.AnnData(X=scipy.sparse.csr_matrix(df, dtype=np.float32),
                               var=pd.DataFrame(index=df.columns),
                               obs=pd.DataFrame(index=df.index.tolist()))
            adata.obs["cellline"] = "8858-1"
    adata.obs["organoid_age_days"] = sample_fn.split("_")[2][1:]
    return adata
