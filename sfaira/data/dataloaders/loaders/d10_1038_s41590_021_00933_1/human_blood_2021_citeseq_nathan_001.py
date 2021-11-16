import anndata
import os
import scipy.sparse
import pandas as pd

disease_map = {"CASE": "pulmonary tuberculosis", "CONTROL": "healthy"}


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read_text(fn, delimiter="\t")
    adata = adata.T
    adata.X = scipy.sparse.csr_matrix(adata.X)
    metadata = pd.read_csv(os.path.join(data_dir, "GSE158769_meta_data.txt.gz"), sep="\t")
    adata.obs = adata.obs.join(metadata.set_index("cell_id"))
    adata.obs["cluster_name"] = adata.obs["cluster_name"].astype("str")

    adata.obs["disease"] = adata.obs["TB_status"].map(disease_map)

    adata.obsm["X_umap"] = adata.obs.loc[:, ["UMAP_1", "UMAP_2"]].values
    return adata
