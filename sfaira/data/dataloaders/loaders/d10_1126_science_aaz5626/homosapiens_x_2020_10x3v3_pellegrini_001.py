import anndata
import os
import scipy.sparse
import pandas as pd
import numpy as np


def load(data_dir, **kwargs):

    em = pd.read_csv(os.path.join(data_dir, "exprMatrix.tsv.gz"), sep="\t", index_col=0).T
    umap_coords = pd.read_csv(os.path.join(data_dir, "Seurat_umap.coords.tsv.gz"), sep="\t", index_col=0, header=None)

    X = scipy.sparse.csr_matrix(em.values, dtype=np.float32)
    var = pd.DataFrame(data={}, index=[i.split("|")[0] for i in em.columns])
    obs = pd.read_csv(os.path.join(data_dir, "meta.tsv"), sep="\t", index_col=0)
    obsm = {"X_umap": umap_coords.loc[obs.index].values.astype(np.float32)}

    assay_diff_dict = {
        "Tel": "Pellegrini, 2020 (doi: 10.1126/science.aaz5626); hCO",
        "ChP1": "Pellegrini, 2020 (doi: 10.1126/science.aaz5626); hChPO",
        "ChP2": "Pellegrini, 2020 (doi: 10.1126/science.aaz5626); hChPO",
        "ChP3": "Pellegrini, 2020 (doi: 10.1126/science.aaz5626); hChPO"
    }

    assay_type_diff_dict = {
        "Tel": "guided",
        "ChP1": "guided",
        "ChP2": "guided",
        "ChP3": "guided"}

    organoid_age_dict = {
        "Tel": "55",
        "ChP1": "27",
        "ChP2": "46",
        "ChP3": "53"}

    obs["assay_diff"] = [assay_diff_dict[i] for i in obs["orig.ident"]]
    obs["assay_type_diff"] = [assay_type_diff_dict[i] for i in obs["orig.ident"]]
    obs["organoid_age_days"] = [organoid_age_dict[i] for i in obs["orig.ident"]]

    adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm)

    return adata
