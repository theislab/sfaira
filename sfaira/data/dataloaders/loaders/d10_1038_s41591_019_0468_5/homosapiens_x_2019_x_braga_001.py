import anndata
import os
import numpy as np
import pandas as pd


def load(data_dir, sample_fn, **kwargs):
    if sample_fn == 'GSE130148':
        fn = [
            os.path.join(data_dir, "GSE130148_raw_counts.csv.gz"),
            os.path.join(data_dir, "GSE130148_barcodes_cell_types.txt.gz"),
        ]
        adata = anndata.read_csv(fn[0]).T
        adata.obs = pd.read_csv(fn[1], sep="\t", index_col=0)

        return adata
    else:
        fn = os.path.join(data_dir, sample_fn)
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)

        return adata
