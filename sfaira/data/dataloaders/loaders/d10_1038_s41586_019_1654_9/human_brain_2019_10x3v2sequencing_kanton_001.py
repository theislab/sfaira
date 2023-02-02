import anndata
import os
import scipy.io
import pandas
import numpy as np


def load(data_dir, **kwargs):
    day_dict = {
        "Organoid-2M": "60",
        "Organoid-4M": "120",
        "Organoid-1M": "30",
        "iPSCs": "0",
        "EB": "4",
        "Neuroectoderm": "10",
        "Neuroepithelium": "15",
    }
    x = scipy.io.mmread(os.path.join(data_dir, 'human_cell_counts_GRCh38.mtx')).T.tocsr().astype(np.float32)
    var = pandas.read_csv(os.path.join(data_dir, 'genes_GRCh38.txt'), sep="\t", index_col=1, names=['ensembl', 'genetype'])
    obs = pandas.read_csv(os.path.join(data_dir, 'metadata_human_cells.tsv'), sep="\t", index_col=0)
    obs["celltype"] = obs['cl_FullLineage'].fillna(obs['cl_LineComp']).fillna("unknown")
    obs["organoid_age_days"] = [day_dict[i] for i in obs["Stage"]]
    adata = anndata.AnnData(X=x, var=var, obs=obs)
    adata.var_names_make_unique()
    return adata
