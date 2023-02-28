import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse
import os
import scanpy as sc


def load(data_dir, **kwargs):
    dic = {
        "sample": ["donor", "organoid_age_days", "protocol"],
        "H1SWeek3": ["H1", "21", "Least Direct [L]"],  # male
        "H1XWeek3": ["H1", "21", "Most Directed [M]"],
        "H1XWeek5": ["H1", "35", "Most Directed [M]"],
        "H1SWeek5": ["H1", "35", "Least Direct [L]"],
        "H1SWeek8": ["H1", "56", "Least Direct [L]"],
        "H1XWeek8": ["H1", "56", "Most Directed [M]"],
        "H1XWeek10": ["H1", "70", "Most Directed [M]"],
        "H1SWeek10": ["H1", "70", "Least Direct [L]"],
        "H1SWeek24": ["H1", "168", "Protocol [L]"],  # male, GEO name is wk24H1S_S2
        "H28126SWeek3": ["H28126", "21", "Least Direct [L]"],  # male
        "H28126XWeek3": ["H28126", "21", "Most Directed [M]"],
        "H28126SWeek5": ["H28126", "35", "Least Direct [L]"],
        "H28126XWeek5": ["H28126", "35", "Most Directed [M]"],
        "H28126PWeek5": ["H28126", "35", "Directed [D]"],
        "H28126S2Week5": ["H28126", "35", "Least Direct [L]"],  # hES-derived cerebral organoid
        "H28126SWeek8": ["H28126", "56", "Least Direct [L]"],
        "H28126XWeek8": ["H28126", "56", "Most Directed [M]"],
        "H28126S2Week8": ["H28126", "56", "Least Direct [L]"],
        "H28126PWeek8": ["H28126", "56", "Directed [D]"],
        "H28126SWeek10": ["H28126", "70", "Least Direct [L]"],
        "H28126XWeek10": ["H28126", "70", "Most Directed [M]"],
        "H28126S2Week10": ["H28126", "70", " Least Direct [L]"],
        "YH10SWeek5": ["WTC10", "35", "Least Direct [L]"],  # male
        "YH10PWeek5": ["WTC10", "35", "Directed [D]"],
        "YH10SWeek8": ["WTC10", "56", "Least Direct [L]"],
        "YH10PWeek8": ["WTC10", "56", "Directed [D]"],
        "YH10PWeek10": ["WTC10", "70", "Directed [D]"],
        "YH10SWeek10": ["WTC10", "70", "Least Direct [L]"],
        "Week5S_1323_4": ["1323_4", "35", "Least Direct [L]"],  # female
        "Week5P_1323_4": ["1323_4", "35", "Directed [D]"],
        "Week8S_1323_4": ["1323_4", "56", "Least Direct [L]"],
        "Week8P_1323_4": ["1323_4", "56", "Directed [D]"],
        "Week10S_1323_4": ["1323_4", "70", "Least Direct [L]"],
        "Week10P_1323_4": ["1323_4", "70", "Directed [D]"],
        "L13234PWeek24": ["1323_4", "168", "Directed [D]"],  # GEO wk2413234P, female; I think this is 1323_4
        "L13234SWeek24": ["1323_4", "168", "Least Direct [L]"],  # GEO wk2413234S, female; I think this is 1323_4
        "L13234SWeek15": ["1323_4", "105", "Least Direct [L]"]  # male, 10x v2, in GEO says its age is 24 weeks, but the name suggests 15 weeks?    
    }
    fn = os.path.join(data_dir, "GSE132672_allorganoids_withnew_matrix.txt.gz")
    df = pd.read_csv(fn, delimiter='\t', compression='gzip', index_col=0)
    adata = ad.AnnData(df.T)
    adata.obs['sample'] = df.columns.str.slice(stop=-17)
    adata.obs[["cell_line", "organoid_age_day", "protocol"]] = [dic[i] for i in adata.obs['sample']]

    return adata
