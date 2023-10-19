import os
import scipy.sparse
import anndata as ad
import numpy as np
import tables


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "pbmc_10k_v3_filtered_feature_bc_matrix.h5")
    with tables.open_file(str(fn), 'r') as f:
        dsets = {}
        for node in f.walk_nodes('/matrix', 'Array'):
            dsets[node.name] = node.read()

    M, N = dsets['shape']
    data = dsets['data']
    if dsets['data'].dtype == np.dtype('int32'):
        data = dsets['data'].view('float32')
        data[:] = dsets['data']
    matrix = scipy.sparse.csr_matrix(
        (data, dsets['indices'], dsets['indptr']),
        shape=(N, M),
    )
    adata = ad.AnnData(
        matrix,
        dict(obs_names=dsets['barcodes'].astype(str)),
        dict(
            var_names=dsets['name'].astype(str),
            gene_ids=dsets['id'].astype(str),
            feature_types=dsets['feature_type'].astype(str),
            genome=dsets['genome'].astype(str),
        ),
    )

    return adata
