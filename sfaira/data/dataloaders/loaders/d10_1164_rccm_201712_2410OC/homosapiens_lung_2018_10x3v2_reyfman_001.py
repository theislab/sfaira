import tarfile

import anndata
import os

import pandas as pd
import scipy.sparse
import h5py


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, 'GSE122960_RAW.tar')
    with tarfile.open(fn) as tar:
        f = h5py.File(tar.extractfile(f'{sample_fn}_filtered_gene_bc_matrices_h5.h5'), 'r')['GRCh38']
        x = scipy.sparse.csc_matrix((f['data'], f['indices'], f['indptr']), shape=f['shape']).T
        var = pd.DataFrame({'feature_id': f['genes'], 'feature_symbol': f['gene_names']})
        obs = pd.DataFrame({'barcode': f['barcodes']})
        adata = anndata.AnnData(X=x, obs=obs, var=var)

    return adata
