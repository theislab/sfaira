import anndata
import tarfile
import gzip
import scipy.io
import os
import pandas as pd


def load(data_dir, sample_fn, **kwargs):
    with tarfile.open(os.path.join(data_dir, 'GSE117770_RAW.tar')) as tar:
        for member in tar.getmembers():
            if "_matrix.mtx.gz" in member.name and sample_fn in member.name:
                name = "_".join(member.name.split("_")[:-1])
                with gzip.open(tar.extractfile(member), "rb") as mm:
                    x = scipy.io.mmread(mm).T.tocsr()
                obs = pd.read_csv(tar.extractfile(name + "_barcodes.tsv.gz"), compression="gzip", header=None,
                                  sep="\t", index_col=0)
                obs.index.name = None
                var = pd.read_csv(tar.extractfile(name + "_genes.tsv.gz"), compression="gzip", header=None,
                                  sep="\t")
                var.columns = ["ensembl", "names"]
                var.index = var["ensembl"].values
                adata = anndata.AnnData(X=x, obs=obs, var=var)
    adata.var_names_make_unique()
    celltypes = pd.read_csv(os.path.join(data_dir, sample_fn + "_annotation.csv"), index_col=0)
    adata = adata[celltypes.index]
    adata.obs["celltypes"] = celltypes

    return adata
