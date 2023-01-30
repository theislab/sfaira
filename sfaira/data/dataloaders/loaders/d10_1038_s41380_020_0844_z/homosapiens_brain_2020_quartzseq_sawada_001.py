import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os
import scanpy as sc


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE120190_RAW.tar")
    adatas = []
    disease = []
    sample = []

    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), compression="gzip", delimiter='\t', header=0, index_col=0)
            d = d.T
            for i in range(d.shape[0]):
                disease.append(member.name.split("_")[2])
                sample.append(member.name.split("_")[2] + "_" + member.name.split("_")[3])
            temp = ad.AnnData(X=scipy.sparse.csr_matrix(d),
                               var=pd.DataFrame(index=d.columns.values),
                               obs=pd.DataFrame(index=d.index))
            adatas.append(temp)

    adata = ad.concat(adatas, join='outer')
    adata.obs['disease'] = disease
    adata.obs['person'] = "DT1_A"
    adata.obs.loc[adata.obs['disease']=="Unaffect", 'person'] = "DT1_U"
    adata.obs['sample'] = sample
    adata.obs['organoid_age_data'] = "35"

    return adata
