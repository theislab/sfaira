import anndata
import gzip
import os
import pandas as pd
import scipy.io
import tarfile


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE164690_RAW.tar")
    with tarfile.open(fn) as tar:
        with gzip.open(tar.extractfile(f"{sample_fn}_matrix.mtx.gz"), "rb") as mm:
            x = scipy.io.mmread(mm).T.tocsr()
        obs = pd.read_csv(tar.extractfile(f"{sample_fn}_barcodes.tsv.gz"), compression="gzip", header=None,
                          sep="\t", index_col=0)
        obs.index.name = None
        var = pd.read_csv(tar.extractfile(f"{sample_fn}_features.tsv.gz"), compression="gzip", header=None, sep="\t")
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
