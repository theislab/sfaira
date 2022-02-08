import anndata
import gzip
import os


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    #with gzip.open(fn, "r") as f:
    adata = anndata.read_h5ad(fn)
    adata.var["feature_types"] = [
        {"ATAC": "peak", "GEX": "rna", "ADT": "protein"}[x]
        for x in adata.var["feature_types"].values
    ]
    return adata
