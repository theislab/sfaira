import os
from typing import Union
import scipy.sparse
import anndata as ad
import numpy as np
import tables

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_blood_2019_10x_10xGenomics_001_unknown"

        self.download_url_data = \
            "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
        self.download_url_meta = None

        self.author = "10x Genomics"
        self.doi = "no_doi_10x_genomics"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "blood"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"
        self.var_ensembl_col = "gene_ids"

    def _load(self):
        fn = os.path.join(self.data_dir, "pbmc_10k_v3_filtered_feature_bc_matrix.h5")
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
            self.adata = ad.AnnData(
                matrix,
                dict(obs_names=dsets['barcodes'].astype(str)),
                dict(
                    var_names=dsets['name'].astype(str),
                    gene_ids=dsets['id'].astype(str),
                    feature_types=dsets['feature_type'].astype(str),
                    genome=dsets['genome'].astype(str),
                ),
            )
