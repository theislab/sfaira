import anndata
import os
import scipy.sparse
import pandas as pd


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE144825_single_CS_GRCh38.txt.gz")
    fn_meta = os.path.join(data_dir, "GSE144825_single_CS_meta.txt.gz")
    data = pd.read_csv(fn, sep="\t", index_col=0).T
    meta = pd.read_csv(fn_meta, sep="\t", index_col=0)
    meta.index.name = None
    meta["Plate"] = meta.index.str.split("_").str[0]
    meta["Plate"] = meta["Plate"].astype("category")
    meta["Line"] = meta["Line"].astype("category")

    adata = anndata.AnnData(
        X=scipy.sparse.csr_matrix(data[data.index.str.endswith("human")].values.copy()),
        var=pd.DataFrame(index=data.columns),
        obs=meta,
        layers={
            "human_counts": scipy.sparse.csr_matrix(data[data.index.str.endswith("human")].values.copy()),
            "chimp_counts": scipy.sparse.csr_matrix(data[data.index.str.endswith("chimp")].values.copy()),
            "total_counts": scipy.sparse.csr_matrix(data[data.index.str.endswith("total")].values.copy()),
        }
    )

    return adata
