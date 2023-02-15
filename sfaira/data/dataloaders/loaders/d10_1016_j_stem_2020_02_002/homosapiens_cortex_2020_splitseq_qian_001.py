import anndata
import os
import tarfile
import pandas as pd
import numpy as np
import scipy.sparse


def load(data_dir, **kwargs):
    with tarfile.open(os.path.join(data_dir, "GSE137941_RAW.tar"), "r") as tar:
        df = pd.read_csv(tar.extractfile(tar.getmembers()[0]), compression='gzip', sep="\t", index_col=0, header=0)
    adata = anndata.AnnData(X=scipy.sparse.csr_matrix(df.T.values, dtype=np.float32),
                            obs=pd.DataFrame(index=df.columns.tolist()),
                            var=pd.DataFrame(index=df.index.tolist()))
    return adata
