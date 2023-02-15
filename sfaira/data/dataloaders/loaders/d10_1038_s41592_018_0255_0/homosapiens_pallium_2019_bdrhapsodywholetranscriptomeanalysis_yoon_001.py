import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os
import scanpy as sc


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE107771_RAW.tar")
    dfs = []
    obs = []
    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), compression="gzip", header=0, index_col=0)
            for i in range(d.shape[0]):
                temp = member.name.split("_")[1] + '_' + member.name.split("_")[2] + member.name.split("_")[3]
                obs.append(temp)
            dfs.append(d)

    adata = ad.AnnData(
        X=scipy.sparse.csr_matrix(np.vstack((dfs[0], dfs[1], dfs[2]))),
        obs=pd.DataFrame({"sample": obs}),
        var=pd.DataFrame(index=dfs[0].columns.values))
    adata.obs['organoid_age_days'] = np.repeat(105, len(adata.obs))

    return adata
