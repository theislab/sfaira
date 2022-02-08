import anndata
import os
import scipy.io
import tarfile
import gzip
import pandas as pd


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE143317_RAW.tar")

    with tarfile.open(fn) as archive:
        file = f'{sample_fn}_filtered_feature_bc_matrix.tar.gz'
        with tarfile.open(fileobj=archive.extractfile(file)) as sample:
            name = sample.getnames()[0]
            with gzip.open(sample.extractfile(name + '/matrix.mtx.gz'), "rb") as mm:
                x = scipy.io.mmread(mm).T.tocsr()
            obs = pd.read_csv(sample.extractfile(name + "/barcodes.tsv.gz"), compression="gzip",
                              header=None, sep="\t", index_col=0)
            obs.index.name = None
            var = pd.read_csv(sample.extractfile(name + "/features.tsv.gz"), compression="gzip",
                              header=None, sep="\t")
            var.columns = ["ensembl", "names", "feature_types"]
            var.index = var["ensembl"].values
            adata = anndata.AnnData(X=x, obs=obs, var=var)
    return adata
