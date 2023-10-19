import numpy as np
import pandas as pd
import anndata as ad
import os
import scipy.io


def load(data_dir, sample_fn, **kwargs):
    sample_fn_dict = {
        "d15": ['counts_org_d15.mtx.gz', 'features_org_d15.tsv.gz', 'meta_org_d15.tsv.gz'],
        "d32to62": ['counts_org1-12.mtx.gz', 'features_org1-12.tsv.gz', 'meta_org1-12.tsv.gz'],
        "microdissected": ['counts.mtx.gz', 'features.tsv.gz', 'meta.tsv.gz'],
    }

    x = scipy.io.mmread(os.path.join(data_dir, sample_fn_dict[sample_fn][0])).T.tocsr().astype(np.float32)
    var = pd.read_csv(os.path.join(data_dir, sample_fn_dict[sample_fn][1]), sep="\t", header=None, index_col=0, names=[None])
    meta = pd.read_csv(os.path.join(data_dir, sample_fn_dict[sample_fn][2]), sep="\t", index_col=0)
    adata = ad.AnnData(X=x, obs=meta, var=var)

    if sample_fn == "d15":
        replace_str = ["is_PSC", "proj_forebrain", "proj_midbrain", "proj_hindbrain", "proj_region", "final_ident"]
        adata.obs[replace_str] = adata.obs[replace_str].astype(str)
        adata.obs["age"] = "15"
        adata.obs["organoid_age_days"] = "15"
    elif sample_fn == "d32to62":
        adata.obs["age"] = [i.split("-d")[1] for i in meta["organoid"]]
        adata.obs["organoid_age_days"] = [i.split("-d")[1] for i in meta["organoid"]]
        adata.obs["celltype"] = pd.read_csv(os.path.join(data_dir, 'cell_annotation_1_12.csv'), index_col=0)["annotation"]
    elif sample_fn == "microdissected":
        adata.obs["age"] = "60"
        adata.obs["organoid_age_days"] = "60"
        adata.obs["celltype"] = pd.read_csv(os.path.join(data_dir, 'cell_annotation_dissected.csv'), index_col=0)["annotation"]

    return adata
