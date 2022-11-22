import anndata
import os
import pandas as pd
import scipy.sparse

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_tar = os.path.join(data_dir, "GSE202297_RAW.tar")
    fn_meta = os.path.join(data_dir, "GSE202297_metadata_q.txt.gz")
    fn_tar = buffered_decompress(fn_tar)
    fn = os.path.join(fn_tar, sample_fn)
    adata = anndata.read_text(fn).transpose()
    adata.X = scipy.sparse.csr_matrix(adata.X)
    tab_meta = pd.read_csv(fn_meta, sep=" ", index_col=False, header=0, compression="gzip")
    del tab_meta["Well_ID"]  # this column de-duplicates otherwise duplicated entries.
    tab_meta = tab_meta.drop_duplicates()
    tab_meta.index = tab_meta["well_coordinates"].values
    adata = adata[[x in tab_meta.index for x in adata.obs_names], :].copy()
    adata.obs = tab_meta.loc[adata.obs_names, :]

    return adata
