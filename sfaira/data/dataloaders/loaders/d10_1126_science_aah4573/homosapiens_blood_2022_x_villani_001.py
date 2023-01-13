import anndata
import os
import pandas as pd


def load(data_dir, **kwargs):
    fn_counts = os.path.join(data_dir, "raw_expression_matrix.txt")
    fn_meta = os.path.join(data_dir, "metadata.txt")
    tab_meta = pd.read_csv(fn_meta, sep="\t", index_col=0, header=0)
    tab_meta = tab_meta.iloc[1:, :].copy()
    adata = anndata.read_text(fn_counts, delimiter="\t").transpose()
    adata.obs_names = [x.replace("_rsem", "") for x in adata.obs_names]
    # Not all cells are annotated:
    adata = adata[[x in tab_meta.index for x in adata.obs_names], :].copy()
    adata.obs = tab_meta
    return adata
