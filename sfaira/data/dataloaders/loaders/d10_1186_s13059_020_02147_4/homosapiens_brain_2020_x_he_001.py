import anndata as ad
import os
import scipy.io
import pandas as pd
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    if sample_fn == "fixation":
        x = scipy.io.mmread(os.path.join(data_dir, 'raw_counts.mtx.gz')).T.tocsr().astype(np.float32)
        var = pd.read_csv(os.path.join(data_dir, 'genes.tsv.gz'), sep="\t", header=None, index_col=0, names=[None])
        obs = pd.read_csv(os.path.join(data_dir, 'cell_metadata.tsv.gz'), sep="\t")
        obs["organoid_age_days"] = "116"
    else:
        x = scipy.io.mmread(os.path.join(data_dir, 'matrix.mtx.gz')).T.tocsr().astype(np.float32)
        var = pd.read_csv(os.path.join(data_dir, 'features.tsv.gz'), sep="\t", header=None, index_col=0, names=[None])
        obs = pd.read_csv(os.path.join(data_dir, 'cells.tsv.gz'), sep="\t")

        meta = pd.read_csv(os.path.join(data_dir, 'meta.tsv.gz'), sep="\t", index_col=0)
        obs["annot_ct"] = meta.loc[obs.index]["annot_ct"]
        obs["organoid_age_days"] = "150"

    adata = ad.AnnData(X=x, obs=obs, var=var)

    return adata
