import anndata
import os


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = adata.layers["counts"]
    adata.obs["donor"] = ["d" + x.split("d")[1] for x in adata.obs["batch"].values]
    adata.obs["site"] = [x.split("d")[0] for x in adata.obs["batch"].values]

    return adata
