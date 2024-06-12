import anndata
import os
import scipy.sparse
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    
    counts = pd.read_csv(data_dir + "GSE115746_cells_exon_counts.csv.gz", index_col=0)
    counts_t = counts.T.copy()

    metadata = pd.read_csv(data_dir + "GSE115746_complete_metadata_28706-cells.csv.gz")
    
    metadata = metadata[metadata.sample_name.isin(counts_t.index.values)].copy()
    metadata.index = metadata.sample_name
    metadata = metadata.reindex(counts_t.index)
    
    
    np.testing.assert_array_equal(metadata.sample_name.values, counts_t.index.values)
    
    adata = AnnData(scipy.sparse.csc_matrix(counts_t.to_numpy()),obs=metadata,)
    adata.var.index = counts_t.columns

    return adata
