import anndata
import os
import pandas as pd
import scipy.io

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, **kwargs):
    fn_barcodes = "GSE136831_AllCells.cellBarcodes.txt.gz"
    fn_counts = "GSE136831_RawCounts_Sparse.mtx.gz"
    fn_features = "GSE136831_AllCells.GeneIDs.txt.gz"
    fn_meta = "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
    fn_counts = buffered_decompress(os.path.join(data_dir, fn_counts))

    tab_barcodes = pd.read_csv(os.path.join(data_dir, fn_barcodes), sep="\t", index_col=0, header=None)
    tab_features = pd.read_csv(os.path.join(data_dir, fn_features), sep="\t", index_col=0)
    tab_meta = pd.read_csv(os.path.join(data_dir, fn_meta), sep="\t", index_col=0)
    x = scipy.io.mmread(fn_counts)
    x = x.transpose()
    adata = anndata.AnnData(x, obs=tab_barcodes, var=tab_features)
    adata.obs = tab_meta
    return adata
