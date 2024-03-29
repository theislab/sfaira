import anndata
import os
import tarfile
import pandas as pd
import scipy.sparse


SAMPLE_DICT = {
    "GSM3589406_PP001swap.filtered.matrix.txt.gz": ["lung", "Donor 1", "healthy"],
    "GSM3589407_PP002swap.filtered.matrix.txt.gz": ["lung", "Donor 1", "stimulated"],
    "GSM3589408_PP003swap.filtered.matrix.txt.gz": ["bone marrow", "Donor 1", "healthy"],
    "GSM3589409_PP004swap.filtered.matrix.txt.gz": ["bone marrow", "Donor 1", "stimulated"],
    "GSM3589410_PP005swap.filtered.matrix.txt.gz": ["lymph node", "Donor 1", "healthy"],
    "GSM3589411_PP006swap.filtered.matrix.txt.gz": ["lymph node", "Donor 1", "stimulated"],
    "GSM3589412_PP009swap.filtered.matrix.txt.gz": ["lung", "Donor 2", "healthy"],
    "GSM3589413_PP010swap.filtered.matrix.txt.gz": ["lung", "Donor 2", "stimulated"],
    "GSM3589414_PP011swap.filtered.matrix.txt.gz": ["bone marrow", "Donor 2", "healthy"],
    "GSM3589415_PP012swap.filtered.matrix.txt.gz": ["bone marrow", "Donor 2", "stimulated"],
    "GSM3589416_PP013swap.filtered.matrix.txt.gz": ["lymph node", "Donor 2", "healthy"],
    "GSM3589417_PP014swap.filtered.matrix.txt.gz": ["lymph node", "Donor 2", "stimulated"],
    "GSM3589418_PP017swap.filtered.matrix.txt.gz": ["blood", "Donor A", "stimulated"],
    "GSM3589419_PP018swap.filtered.matrix.txt.gz": ["blood", "Donor A", "healthy"],
    "GSM3589420_PP019swap.filtered.matrix.txt.gz": ["blood", "Donor B", "stimulated"],
    "GSM3589421_PP020swap.filtered.matrix.txt.gz": ["blood", "Donor B", "healthy"],
}


def load(data_dir, sample_fn, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE126030_RAW.tar"),
        os.path.join(data_dir, "donor1.annotation.txt"),
        os.path.join(data_dir, "donor2.annotation.txt")
    ]
    with tarfile.open(fn[0]) as tar:
        df = pd.read_csv(tar.extractfile(sample_fn), compression="gzip", sep="\t")
    df.index = [i.split(".")[0] for i in df["Accession"]]
    var = pd.concat([df.pop(x) for x in ["Gene", "Accession"]], 1)
    if df.columns[-1].startswith("Un"):
        df.drop(df.columns[len(df.columns) - 1], axis=1, inplace=True)
    adata = anndata.AnnData(df.T)
    adata.var = var
    adata.obs["donor"] = SAMPLE_DICT[sample_fn][1]
    adata.obs.index = sample_fn.split("_")[1].split("s")[0] + "nskept." + adata.obs.index
    adata.obs["cell_ontology_class"] = "unknown"
    df1 = pd.read_csv(fn[1], sep="\t", index_col=0, header=None)
    df2 = pd.read_csv(fn[2], sep="\t", index_col=0, header=None)
    for i in df1.index:
        adata.obs["cell_ontology_class"].loc[i] = df1.loc[i][1]
    for i in df2.index:
        adata.obs["cell_ontology_class"].loc[i] = df2.loc[i][1]
    adata.X = scipy.sparse.csc_matrix(adata.X)

    return adata
