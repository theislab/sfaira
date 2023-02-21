import anndata as ad
import os
import numpy as np
import pandas as pd
import scipy.sparse
import gzip
import io


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE115011_marton_all_cells.csv.gz")
    d = pd.read_csv(fn, index_col=0).T
    adata = ad.AnnData(X=scipy.sparse.csr_matrix(d).astype(np.float32),
                       var=pd.DataFrame(index=d.columns.values),
                       obs=pd.DataFrame(index=d.index.str.strip()))
    df_str = ""
    with gzip.open(os.path.join(data_dir, "GSE115011_series_matrix.txt.gz"), "rb") as m:
        for line in m.readlines():
            linestr = line.decode('utf8')
            if linestr.startswith("!Sample"):
                df_str += linestr
    df = pd.read_csv(io.StringIO(df_str), sep="\t", header=0, index_col=0).T
    df.rename(index={'Rebecca-Oligo_HMMKJ_L7_AAGAGGCA-ACTCTAGG': 'SCGPM_Rebecca-Oligo_HMMKJ_L7_AAGAGGCA-ACTCTAGG'}, inplace=True)
    adata.obs["cellline"] = pd.Categorical(df.loc[adata.obs.index].iloc[:, 10].str.split().str[1].values)
    adata.obs["organoid_age_days"] = "127"
    adata.var_names_make_unique()

    return adata
