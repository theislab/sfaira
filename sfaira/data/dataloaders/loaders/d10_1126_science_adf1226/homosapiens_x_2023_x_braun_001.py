import os
import numpy as np
import anndata


def load(data_dir, **kwargs):
    age_dict = {
        'BRC2006': '8.0',
        'BRC2021': '8.0',
        'BRC2057': '8.1',
        'BRC2061': '6.9',
        'BRC2073': '6.6',
        'BRC2106A': '6.6',
        'BRC2110': '6.9',
        'BRC2114': '6.9',
        'BRC2147': '6.7',
        'BRC2191': '6.9',
        'XDD:313': '8.5',
        'XDD:326': '6.0',
        'XDD:334': '8.0',
        'XDD:342': '8.5',
        'XDD:348': '5.0',
        'XDD:351': '12.0',
        'XDD:358': '11.5',
        'XDD:359': '13.0',
        'XDD:385': '14.0',
        'XDD:395': '6.0',
        'XDD:398': '7.0',
        'XDD:400': '5.5',
        'XHU:292': '9.5',
        'XHU:297': '10.0',
        'XHU:305': '7.5',
        'XHU:307': '9.2'
    }

    adata = anndata.read_h5ad(os.path.join(data_dir, 'human_dev_GRCh38-3.0.0.h5ad'))

    adata.obs.index = adata.obs.index.str[2:-1]
    for c in adata.obs.loc[:, adata.obs.dtypes == "category"].columns:
        if all(adata.obs[c].cat.categories.str.startswith("b'")):
            adata.obs[c] = adata.obs[c].astype(str).str[2:-1].astype("category")
    adata.var_names_make_unique()

    adata.obs["Age"] = adata.obs["donor_id"].replace(age_dict)
    adata.obs["age_days"] = np.rint(np.multiply(np.array(adata.obs["Age"].astype(np.float32)), 7.))

    return adata
