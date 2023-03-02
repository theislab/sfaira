import numpy as np
import pandas as pd
import anndata as ad
import os
import scipy.io


def load(data_dir, **kwargs):
    x0 = scipy.io.mmread(os.path.join(data_dir, 'counts.mtx.gz')).T.tocsr().astype(np.float32)
    var0 = pd.read_csv(os.path.join(data_dir, 'features.tsv.gz'), sep="\t", header=None, index_col=0, names=[None])
    meta0 = pd.read_csv(os.path.join(data_dir, 'meta.tsv.gz'), sep="\t", index_col=0)
    meta0['organoid_age_days'] = "60"
    meta0['experiment'] = "E-MTAB-10973"
    meta0['sample'] = meta0['orig.ident']
    adata0 = ad.AnnData(X=x0, obs=meta0, var=var0)

    x1 = scipy.io.mmread(os.path.join(data_dir, 'counts_org1-12.mtx.gz')).T.tocsr().astype(np.float32)
    var1 = pd.read_csv(os.path.join(data_dir, 'features_org1-12.tsv.gz'), sep="\t", header=None, index_col=0, names=[None])
    meta1 = pd.read_csv(os.path.join(data_dir, 'meta_org1-12.tsv.gz'), sep="\t", index_col=0)
    meta1['organoid_age_days'] = [i.split("d")[1] for i in meta1["organoid"]]
    meta1['sample'] = meta1['organoid']
    meta1['experiment'] = "E-MTAB-10974_org1-12"
    adata1 = ad.AnnData(X=x1, obs=meta1, var=var1)

    x2 = scipy.io.mmread(os.path.join(data_dir, 'counts_org_d15.mtx.gz')).T.tocsr().astype(np.float32)
    var2 = pd.read_csv(os.path.join(data_dir, 'features_org_d15.tsv.gz'), sep="\t", header=None, index_col=0, names=[None])
    meta2 = pd.read_csv(os.path.join(data_dir, 'meta_org_d15.tsv.gz'), sep="\t", index_col=0)
    meta2['experiment'] = "E-MTAB-10974_org_d15"
    adata2 = ad.AnnData(X=x2, obs=meta2, var=var2)
    adata2.obs['organoid_age_days'] = "15"

    adata = adata0.concatenate([adata1, adata2])

    return adata
