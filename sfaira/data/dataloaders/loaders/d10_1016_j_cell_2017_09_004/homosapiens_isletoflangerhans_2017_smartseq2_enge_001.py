import os
import tarfile
import gzip
from io import StringIO
import anndata as ad
import pandas as pd
import scipy.sparse


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE81547_RAW.tar"),
        os.path.join(data_dir, "GSE81547_series_matrix.txt.gz")
    ]
    dfs = []
    with tarfile.open(fn[0]) as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), compression="gzip", header=None, sep="\t", index_col=0,
                            names=[member.name.split("_")[0]])
            dfs.append(d)
    adata = ad.AnnData(pd.concat(dfs, axis=1).iloc[1:-6].T)
    adata.X = scipy.sparse.csc_matrix(adata.X)
    with gzip.open(fn[1]) as f:
        file_content = [i.decode("utf-8") for i in f.readlines()]
    inputstring = ""
    for line in file_content:
        if "ID_REF" in line:
            inputstring += line
        if "!Sample_title" in line:
            inputstring += line[1:]
        if "!Sample_characteristics_ch1\t\"inferred_cell_type: alpha" in line:
            inputstring += line[1:]
    data = StringIO(inputstring)
    d = pd.read_csv(data, sep="\t").T
    d.columns = d.iloc[0]
    d.drop("Sample_title", inplace=True)
    d = d.reset_index().set_index("ID_REF")
    d.columns.name = None
    d.index.name = None
    adata.obs["celltype"] = [d.loc[i]["Sample_characteristics_ch1"].split(": ")[1] for i in adata.obs.index]
    adata.obs["patient"] = ["_".join(d.loc[i]["index"].split("_")[:2]) for i in adata.obs.index]

    return adata
