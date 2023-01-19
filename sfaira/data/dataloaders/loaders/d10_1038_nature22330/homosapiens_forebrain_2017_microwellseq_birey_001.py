import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os
# import any packages you require for dataloading here. you can assume packages like numpy and pandas being available


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir, **kwargs):
    # replace my-data-file.h5ad with the filename you're loading
    fn = os.path.join(data_dir, "GSE93811_RAW.tar")
    # replace the simple data loading code below with the code required to load your data file(s)  
    dfs = []
    obs = []
    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), compression="gzip", header=0, index_col=None)
            for i in range(d.shape[0]):
                obs.append(member.name.split("_")[1])
            dfs.append(d)
    adata = ad.AnnData(X=scipy.sparse.csr_matrix(np.vstack((dfs[0], dfs[1]))), obs=pd.DataFrame({"sample": obs}), var=pd.DataFrame(index=dfs[0].columns.values))

    return adata  # your load function needs to return an AnnData object
