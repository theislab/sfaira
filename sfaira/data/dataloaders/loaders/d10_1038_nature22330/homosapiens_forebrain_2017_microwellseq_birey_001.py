import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE93811_RAW.tar")
    dfs = []
    obs = []
    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), compression="gzip", header=0, index_col=None)
            for i in range(d.shape[0]):
                obs.append(member.name.split("_")[1])
            dfs.append(d)
    adata = ad.AnnData(X=scipy.sparse.csr_matrix(np.vstack((dfs[0], dfs[1]))), obs=pd.DataFrame({"sample": obs}), var=pd.DataFrame(index=dfs[0].columns.values))
    adata.obs['organoid_age_days'] = "105"
    adata.obs['cell_line'] = np.repeat('2242-1', len(adata.obs))
    return adata
