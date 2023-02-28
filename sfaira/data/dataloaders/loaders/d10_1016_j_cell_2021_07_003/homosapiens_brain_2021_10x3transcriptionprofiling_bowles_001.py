import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os
import scanpy as sc


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE171345_RAW.tar")
    dfs = []
    obs = []
    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), delimiter='\t', compression="gzip", header=0, index_col=None)
            for i in range(d.shape[0]):
                obs.append(member.name.split("_")[0])
            dfs.append(d)

    dic = {
        "sample" : ["cell_line", "organoid_age_days", "genotype"],
        "GSM5223939": ["NDB06 & NDB09", "2mo", "V337V & V337M"], 
        "GSM5223940": ["NDB06 & NDB09", "2mo", "V337V & V337M"],
        "GSM5223941": ["NDB06", "4mo", "V337V"],
        "GSM5223942": ["NDB09", "4mo", "V337M"],
        "GSM5223943": ["GIH6-E11 & GIH6-A02", "2mo", "V337V & V337M"],
        "GSM5223944": ["GIH7-B12 & GIH7-F02", "2mo", "V337V"],
        "GSM5223945": ["NDB06 & NDB09 & GIH7-A01", "2mo", "V337V & V337M"],
        "GSM5223946": ["NDB06 & NDB09 & GIH7-A01", "2mo", "V337V & V337M"],
        "GSM5223947": ["GIH6-A02 & GIH7-A01", "2mo", "V337M"],
        "GSM5223948": ["GIH7-B12 & GIH7-F02", "2mo", "V337V"],
        "GSM5223949": ["NDB06 & NDB09", "2mo", "V337V & V337M"],
        "GSM5223950": ["GIH6-E11 & GIH7-B12", "2mo", "V337V"],
        "GSM5223951": ["GIH7-F02 & GIH6-A02", "2mo", "V337V & V337M"],
        "GSM5223952": ["GIH7-A01 & NDB06 & NDB09", "2mo", "V337V & V337M"],
        "GSM5223953": ["GIH6-E11", "2mo", "V337V"],
        "GSM5223954": ["NDB06 & NDB09", "6mo", "V337V & V337M"],
        "GSM5223955": ["GIH6-E11 & GIH6-A02", "4mo", "V337V & V337M"],
        "GSM5223956": ["GIH7-B12 & GIH7-F02", "4mo", "V337V"],
        "GSM5223957": ["NDB06 & NDB09", "4mo", "V337V & V337M"],
        "GSM5223958": ["GIH7-A01", "4mo", "V337M"], # mutant tau associated with frontotemporal dementia
        "GSM5223959": ["GIH6-E11 & GIH7-B12", "4mo", "V337V"], # normal tau
        "GSM5223960": ["NDB06 & NDB09", "4mo", "V337V & V337M"],
        "GSM5223961": ["GIH7-F02 & GIH6-A02", "4mo", "V337V & V337M"],
        "GSM5223962": ["GIH7-A01 & GIH6-A02", "4mo", "V337M"],
        "GSM5223963": ["GIH6-E11 & GIH7-A01", "4mo", "V337V & V337M"],
        "GSM5223964": ["GIH7-B12 & GIH7-F02", "4mo", "V337V"],
        "GSM5223965": ["GIH6-E11 & GIH6-A02", "6mo", "V337V & V337M"],
        "GSM5223966": ["GIH7-B12 & GIH7-F02", "6mo", "V337V"],
        "GSM5223967": ["NDB06 & GIH7-A01", "6mo", "V337V & V337M"],
        "GSM5223968": ["NDB06 & NDB09 & GIH7-A01", "6mo", "V337V & V337M"],
        "GSM5223969": ["GIH6-A02 & GIH7-B12", "6mo", "V337V & V337M"],
        "GSM5223970": ["GIH7-F02 & GIH7-B12", "6mo", "V337V"],
        "GSM5223971": ["NDB06 & GIH7-F02", "6mo", "V337V"],
        "GSM5223972": ["GIH7-A01 & GIH6-E11", "6mo", "V337V & V337M"],
        "GSM5223973": ["GIH6-A02 & NDB09", "6mo", "V337M"],
        "GSM5223974": ["GIH7-A01 & GIH6-E11", "6mo", "V337V & V337M"]
    }

    adata = ad.AnnData(
        X=scipy.sparse.csr_matrix(np.vstack(dfs)), 
        obs=pd.DataFrame({"sample": obs}), 
        var=pd.DataFrame(index=dfs[0].columns.values))
    adata.obs[["cell_line", "organoid_age_days", "genotype"]] = [dic[i] for i in adata.obs['sample']]

    fn_meta = os.path.join(data_dir, "meta.tsv.gz")
    meta = pd.read_csv(fn_meta, delimiter='\t', compression='gzip')
    meta = meta.set_index(meta['CellID'].str.split('_', expand=True)[0])

    return adata  