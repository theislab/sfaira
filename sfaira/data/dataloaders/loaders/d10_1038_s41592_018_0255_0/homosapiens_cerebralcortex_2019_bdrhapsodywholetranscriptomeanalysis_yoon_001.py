import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE107771_RAW.tar")
    with tarfile.open(fn, "r") as tar:
        df = pd.read_csv(tar.extractfile(f"{sample_fn}_MolsPerCell.csv.gz"), compression="gzip", index_col=0)
        adata = ad.AnnData(X=scipy.sparse.csr_matrix(df, dtype=np.float32),
                           var=pd.DataFrame(index=df.columns),
                           obs=pd.DataFrame(index=df.index.tolist()))
    adata.obs['organoid_age_days'] = "105"

    return adata
