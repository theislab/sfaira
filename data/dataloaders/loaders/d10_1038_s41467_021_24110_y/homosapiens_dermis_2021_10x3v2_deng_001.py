import anndata
import gzip
import os
import pandas as pd
import scipy.io
import tarfile

sample_dict = {
    'GSM4994379_KL1': "KF1",
    'GSM4994380_KL2': "KF2",
    'GSM4994381_KL3': "KF3",
    'GSM4994382_NS1': "NF1",
    'GSM4994383_NS2': "NF2",
    'GSM4994384_NS3': "NF3",
}

disease_dict = {
    'GSM4994379_KL1': "keloid",
    'GSM4994380_KL2': "keloid",
    'GSM4994381_KL3': "keloid",
    'GSM4994382_NS1': "healthy",
    'GSM4994383_NS2': "healthy",
    'GSM4994384_NS3': "healthy",
}

sex_dict = {
    'GSM4994379_KL1': "male",
    'GSM4994380_KL2': "male",
    'GSM4994381_KL3': "female",
    'GSM4994382_NS1': "male",
    'GSM4994383_NS2': "male",
    'GSM4994384_NS3': "female",
}

dev_dict = {
    'GSM4994379_KL1': "20-year-old human stage",
    'GSM4994380_KL2': "23-year-old human stage",
    'GSM4994381_KL3': "34-year-old human stage",
    'GSM4994382_NS1': "39-year-old human stage",
    'GSM4994383_NS2': "28-year-old human stage",
    'GSM4994384_NS3': "26-year-old human stage",
}


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE163973_RAW.tar")
    with tarfile.open(fn) as tar:
        with tarfile.open(fileobj=tar.extractfile(sample_fn+"_matrix.tar.gz")) as tar2:
            with gzip.open(tar2.extractfile(f"{tar2.getnames()[0]}/matrix.mtx.gz"), "rb") as mm:
                x = scipy.io.mmread(mm).T.tocsr()
            obs = pd.read_csv(tar2.extractfile(f"{tar2.getnames()[0]}/barcodes.tsv.gz"), compression="gzip", header=None, sep="\t", index_col=0)
            obs.index.name = None
            obs.index = [f"{i.split('-')[0]}_{sample_dict[sample_fn]}" for i in obs.index]
            var = pd.read_csv(tar2.extractfile(f"{tar2.getnames()[0]}/features.tsv.gz"), compression="gzip", header=None, sep="\t")
            var.columns = ["ensembl", "symbol", "feature_class"]
            var.index = var["symbol"].values
            adata = anndata.AnnData(X=x, obs=obs, var=var)

            meta = pd.read_csv(os.path.join(data_dir, "GSE163973_integrate.all.NS.all.KL_cell.meta.data.csv.gz", index_col=0)
            meta.index = [f"{ind.split('_')[0]}_{meta['orig.ident'].iloc[j]}" for j, ind in enumerate(meta.index)]
            meta = meta[meta["orig.ident"] == sample_dict[sample_fn]]
            adata = adata[meta.index].copy()
            adata.obs = meta
            adata.obs["disease"] = disease_dict[sample_fn]
            adata.obs["sex"] = sex_dict[sample_fn]
            adata.obs["age"] = dev_dict[sample_fn]

            return adata
