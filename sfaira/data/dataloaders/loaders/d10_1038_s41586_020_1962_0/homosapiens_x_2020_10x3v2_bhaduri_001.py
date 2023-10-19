import os
import pandas as pd
import numpy as np
import scipy.sparse
import gzip
import anndata as ad


def load(data_dir, sample_fn, **kwargs):
    def read_txt_sparse(fn, sep="\t"):
        X_data = []
        X_row, X_col = [], []
        index = []
        with gzip.open(fn, mode="rt") as f:
            header = f.readline().strip().split(sep)
            for row_idx, line in enumerate(f.readlines()):
                row = line.split(sep)
                index.append(row.pop(0))
                row = np.array(row, dtype=np.float32)
                col_inds, = np.nonzero(row)
                X_col.extend(col_inds)
                X_row.extend([row_idx] * len(col_inds))
                X_data.extend(row[col_inds])
        X = scipy.sparse.coo_matrix((X_data, (X_row, X_col)), shape=(len(index), len(header)), dtype=np.float32).T.tocsr()
        return X, header, index

    if sample_fn == "organoids":
        fn = os.path.join(data_dir, 'organoidsmatrix_nonorm.txt.gz')
        X, header, index = read_txt_sparse(fn)

        adata = ad.AnnData(X=X, obs=pd.DataFrame(index=header), var=pd.DataFrame(index=index))
        meta = pd.read_excel(os.path.join(data_dir, "41586_2020_1962_MOESM3_ESM.xlsx"),
                             sheet_name="STable 1 Cell Metadata", usecols="N:AB", skiprows=1, index_col=0)
        meta.columns = [i.strip(".1") for i in meta.columns]

        adata.obs = adata.obs.reset_index()["index"].str.split(" ", expand=True)
        adata.obs = pd.DataFrame(index=adata.obs[1] + "_" + adata.obs[0])
        new_index = list(set(adata.obs.index) & set(meta.index))
        adata = adata[new_index].copy()
        adata.obs = meta.loc[new_index].copy()
        adata.obs["assay_type_diff"] = adata.obs["Protocol"].replace({"Least Directed": "guided", "Directed": "guided", "Most Directed": "guided"})
        adata.obs["assay_diff"] = adata.obs["Protocol"].replace({
            "Least Directed": "Velasco, 2019 (doi: 10.1038/s41586-019-1289-x)",
            "Directed": "Bhaduri, 2020 (doi: 10.1038/s41586-020-1962-0); directed",
            "Most Directed": "Bhaduri, 2020 (doi: 10.1038/s41586-020-1962-0); most directed"
        })
    else:
        fn = os.path.join(data_dir, 'primarymatrix_nonorm.txt.gz')
        X, header, index = read_txt_sparse(fn)

        adata = ad.AnnData(X=X, obs=pd.DataFrame(index=header), var=pd.DataFrame(index=index))
        meta = pd.read_csv(os.path.join(data_dir, 'meta.tsv'), sep="\t", index_col=0)
        tsne = pd.read_csv(os.path.join(data_dir, 'Seurat_tsne.coords.tsv.gz'), sep="\t", index_col=0, header=None)

        adata.obs = adata.obs.reset_index()["index"].str.split(" ", expand=True)
        adata.obs[2] = adata.obs[1].str[:4]
        adata.obs[3] = adata.obs[1].str[4:].str.lower().replace({"somato": "somatosensory"})
        adata.obs.index = adata.obs[[2, 3, 0]].apply(lambda row: '_'.join(row), axis=1).values

        meta["arealower"] = meta["Area"].str.lower().str.split().str[0]
        meta["barcodes"] = meta.index.str.extract('([A-Z]{16})')[0].to_list()
        meta.index = meta[["Individual", "arealower", "barcodes"]].apply(lambda row: '_'.join(row), axis=1).values
        tsne.index = meta[["Individual", "arealower", "barcodes"]].apply(lambda row: '_'.join(row), axis=1).values
        meta = meta.loc[meta.index.drop_duplicates(keep=False)].copy()
        tsne = tsne.loc[tsne.index.drop_duplicates(keep=False)].copy()

        consensus_index = list(set(meta.index) & set(adata.obs.index))
        adata = adata[consensus_index].copy()
        adata.obs = meta.loc[consensus_index].copy()
        adata.obsm["X_tsne"] = tsne.loc[consensus_index].values.astype(np.float32)

    adata.obs["organoid_age_days"] = adata.obs["Age"] * 7
    adata.obs["organoid_age_days"] = adata.obs["organoid_age_days"].astype("str")
    adata.obs["Age"] = adata.obs["Age"].astype("str")
    adata.obs["celltype"] = adata.obs["Type"] + "_" + adata.obs["Subtype"]
    adata.var_names_make_unique()

    return adata
