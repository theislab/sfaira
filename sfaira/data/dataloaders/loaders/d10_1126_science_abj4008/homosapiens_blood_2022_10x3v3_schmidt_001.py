import os

import pandas as pd
import scanpy as sc


def load(data_dir, sample_fn, **kwargs):
    adata = sc.read_10x_mtx(data_dir, prefix="GSE190604_")
    fn_meta = os.path.join(data_dir, "GSE190604_cellranger-guidecalls-aggregated-unfiltered.txt.gz")
    tab_meta = pd.read_csv(fn_meta, compression="gzip", sep="\t")
    tab_meta.index = tab_meta["cell_barcode"].values
    del tab_meta["cell_barcode"]
    adata.obs = pd.concat([adata.obs, tab_meta], axis=1)
    return adata
