import anndata
import os


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = adata.layers["counts"]
    adata.var["feature_types"] = [
        {"GEX": "rna", "ADT": "protein"}[x]
        for x in adata.var["feature_types"].values
    ]
    return adata
