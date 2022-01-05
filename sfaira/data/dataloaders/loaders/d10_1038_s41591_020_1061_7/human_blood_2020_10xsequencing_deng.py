
import anndata
import os
import pandas as pd
import scipy.io
import tarfile

sex_map = {"M": "male", "F": "female"}
disease_map = {"DLBCL": "diffuse large B-cell lymphoma",
               "tFL": "follicular lymphoma",
               "PMBCL": "primary mediastinal large B-cell lymphoma"}


def load(data_dir, sample_fn, **kwargs):
    fn = [os.path.join(data_dir, "GSE151511_RAW.tar"),
          os.path.join(data_dir, "41591_2020_1061_MOESM3_ESM.xlsx")]

    with tarfile.open(fn[0]) as tar:
        with tarfile.open(fileobj=tar.extractfile(sample_fn)) as sample:
            name = sample.getnames()[0]

            x = scipy.io.mmread(sample.extractfile(f'{name}/{name}_matrix.mtx')).T.tocsr()
            obs = pd.read_csv(sample.extractfile(f'{name}/{name}_barcodes.tsv'), header=None, sep="\t", index_col=0)
            obs.index.name = None
            var = pd.read_csv(sample.extractfile(f'{name}/{name}_genes.tsv'), header=None, sep="\t")
            var.columns = ["ensembl", "names"]
            var.index = var["ensembl"].values
            adata = anndata.AnnData(X=x, obs=obs, var=var)
    metadata = pd.read_excel(fn[1], skiprows=2, skipfooter=3)

    adata.obs['Sample ID'] = name
    adata.obs = adata.obs.reset_index().merge(metadata, on='Sample ID').set_index("index")
    adata.obs["Age"] = [f"{x}-year-old human stage" for x in adata.obs["Age"].values]
    adata.obs["Gender"] = adata.obs["Gender"].map(sex_map)
    adata.obs["Histology"] = adata.obs["Histology"].map(disease_map)

    adata.obs = adata.obs.astype(str)
    return adata
