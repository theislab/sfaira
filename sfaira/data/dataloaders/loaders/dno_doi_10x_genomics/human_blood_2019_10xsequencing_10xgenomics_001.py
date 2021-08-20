import os
import scipy.sparse
import anndata as ad
import numpy as np
import tables

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = \
            "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v3"
        self.author = "10x Genomics"
        self.disease = "healthy"
        self.doi_journal = "no_doi_10x_genomics"
        self.normalization = "raw"
        self.organ = "blood"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.year = 2019

        self.gene_id_symbols_var_key = "index"
        self.gene_id_ensembl_var_key = "gene_ids"

        self.set_dataset_id(idx=1)


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
