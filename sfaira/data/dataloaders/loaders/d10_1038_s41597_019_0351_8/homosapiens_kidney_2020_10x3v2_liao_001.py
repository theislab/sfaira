import anndata
import os
import pandas as pd
import scipy.io
import gzip
import tarfile


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE131685_RAW.tar")
    adatas = []
    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            if "_matrix.mtx.gz" in member.name:
                name = "_".join(member.name.split("_")[:-1])
                with gzip.open(tar.extractfile(member), "rb") as mm:
                    X = scipy.io.mmread(mm).T.tocsr()
                obs = pd.read_csv(tar.extractfile(name + "_barcodes.tsv.gz"), compression="gzip", header=None,
                                  sep="\t", index_col=0)
                obs.index.name = None
                var = pd.read_csv(tar.extractfile(name + "_features.tsv.gz"), compression="gzip", header=None,
                                  sep="\t").iloc[:, :2]
                var.columns = ["ensembl", "names"]
                var.index = var["ensembl"].values
                adata = anndata.AnnData(X=X, obs=obs, var=var)
                adata.obs["sample"] = name
                adatas.append(adata)
    adata = adatas[0].concatenate(adatas[1:])

    return adata
