import anndata
import os
from typing import Union

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader requires manual preprocessing of the raw datafile. To download the data, use the link in the
    `.download_website` attribute of this class. To create the file required by this dataloader, run the following
    python code:

    import scanpy
    scanpy.read_10x_h5('pbmc_10k_v3_filtered_feature_bc_matrix.h5').write('pbmc_10k_v3_filtered_feature_bc_matrix.h5ad')

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_blood_2019_10x_10xGenomics_001_unknown"

        self.download = "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
        self.download_meta = None

        self.author = '10x Genomics'
        self.doi = "d_nan"
        self.healthy = True
        self.normalization = 'raw'
        self.organ = "blood"
        self.organism = "human"
        self.protocol = '10x'
        self.state_exact = 'healthy'
        self.year = 2019

        self.var_symbol_col = 'index'
        self.var_ensembl_col = 'gene_ids'

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "blood", "pbmc_10k_v3_filtered_feature_bc_matrix.h5ad")
        self.adata = anndata.read(fn)
