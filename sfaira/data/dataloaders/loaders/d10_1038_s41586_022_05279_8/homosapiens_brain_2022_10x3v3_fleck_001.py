import anndata as ad
import os
import scipy.io
import pandas as pd
import numpy as np


def load(data_dir, **kwargs):
    x = scipy.io.mmread(os.path.join(data_dir, 'counts.mtx.gz')).T.tocsr().astype(np.float32)
    var = pd.read_csv(os.path.join(data_dir, 'features.tsv.gz'), sep="\t", header=None, index_col=0, names=[None])
    obs = pd.read_csv(os.path.join(data_dir, 'barcodes.tsv.gz'), sep="\t", index_col=0, header=None, names=[None])

    meta = pd.read_csv(os.path.join(data_dir, 'meta.tsv.gz'), sep="\t", index_col=0)
    meta = meta.loc[obs.index]
    meta["organoid_age_days"] = [i.split("_")[1] for i in meta["sample"]]
    meta["celltype"] = meta["lineage"] + "-" + meta["stage_manual"]

    adata = ad.AnnData(X=x, obs=meta, var=var)

    return adata
