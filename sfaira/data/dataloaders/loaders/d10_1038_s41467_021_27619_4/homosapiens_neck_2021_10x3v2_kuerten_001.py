import anndata
import gzip
import os
import pandas as pd
import scipy.io

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE164690_RAW.tar")
    tar_dir = buffered_decompress(fn)
    fn_mtx = buffered_decompress(os.path.join(tar_dir, f"{sample_fn}_matrix.mtx.gz"))
    x = scipy.io.mmread(fn_mtx).T.tocsr()
    fn_obs = buffered_decompress(os.path.join(tar_dir, f"{sample_fn}_barcodes.tsv.gz"))
    obs = pd.read_csv(fn_obs, header=None, sep="\t", index_col=0)
    obs.index.name = None
    fn_var = buffered_decompress(os.path.join(tar_dir, f"{sample_fn}_features.tsv.gz"))
    var = pd.read_csv(fn_var, header=None, sep="\t")
    var.columns = ["ensembl", "symbol", "feature_class"]
    var.index = var["ensembl"].values
    adata = anndata.AnnData(X=x, obs=obs, var=var)
    # Annotate organism in .var
    adata.var["organism"] = ["Homo sapiens" if x.startswith("ENSG") else "Human papillomavirus"
                             for x in adata.var["ensembl"].values]
    # TODO: allow HPV back in once multi-organism is supported.
    adata = adata[:, adata.var["organism"].values == "Homo sapiens"].copy()
    adata.obs["patient"] = sample_fn.split("_")[1]
    return adata
