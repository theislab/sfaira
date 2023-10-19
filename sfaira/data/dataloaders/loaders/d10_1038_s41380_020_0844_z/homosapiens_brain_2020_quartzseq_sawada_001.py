import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE120190_RAW.tar")
    with tarfile.open(fn) as tar:
        d = pd.read_csv(tar.extractfile(f"{sample_fn}_63mer_seqlev_d2_dge.txt.gz"), compression="gzip", delimiter='\t', header=0, index_col=0).T
        adata = ad.AnnData(X=scipy.sparse.csr_matrix(d).astype(np.float32),
                           var=pd.DataFrame(index=d.columns.values),
                           obs=pd.DataFrame(index=d.index))
    adata.obs['organoid_age_days'] = "35"

    return adata
