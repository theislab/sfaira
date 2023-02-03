import anndata
import os
import scipy.sparse
import pandas as pd
import numpy as np


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE144825_single_CS_GRCh38.txt.gz")
    fn_meta = os.path.join(data_dir, "GSE144825_single_CS_meta.txt.gz")
    data = pd.read_csv(fn, sep="\t", index_col=0).T
    meta = pd.read_csv(fn_meta, sep="\t", index_col=0)
    meta.index.name = None
    meta["Plate"] = meta.index.str.split("_").str[0]
    meta["organoid_age_days"] = "150"
    meta["organoid_age_days"] = meta["organoid_age_days"].astype("category")
    meta["Plate"] = meta["Plate"].astype("category")
    meta["Line"] = meta["Line"].astype("category")

    adata = anndata.AnnData(
        X=scipy.sparse.csr_matrix(data[data.index.str.endswith("human")].values.copy(), dtype=np.float32),
        var=pd.DataFrame(index=data.columns),
        obs=meta,
        layers={
            "human_counts": scipy.sparse.csr_matrix(data[data.index.str.endswith("human")].values.copy(), dtype=np.float32),
            "chimp_counts": scipy.sparse.csr_matrix(data[data.index.str.endswith("chimp")].values.copy(), dtype=np.float32),
            "total_counts": scipy.sparse.csr_matrix(data[data.index.str.endswith("total")].values.copy(), dtype=np.float32),
        }
    )

    return adata
