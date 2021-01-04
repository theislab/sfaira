import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class.

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
        self.id = "human_eye_2019_10x_lukowski_001_10.15252/embj.2018100811"

        self.download = "https://covid19.cog.sanger.ac.uk/lukowski19.processed.h5ad"
        self.download_meta = None

        self.author = 'Wong'
        self.doi = '10.15252/embj.2018100811'
        self.healthy = True
        self.normalization = 'raw'
        self.organ = "eye"
        self.organism = "human"
        self.protocol = '10x'
        self.state_exact = 'healthy'
        self.sub_tissue = "retina"
        self.year = 2019

        self.var_symbol_col = 'index'
        self.var_ensembl_col = 'gene_ids'

        self.obs_key_cellontology_original = 'CellType'

        self.class_maps = {
            "0": {
                'Muller cell': 'Muller cell',
                'amacrine cell': 'Amacrine cell',
                'microglial cell': 'Microglia',
                'retinal bipolar neuron type A': 'Retinal bipolar neuron type A',
                'retinal bipolar neuron type B': 'Retinal bipolar neuron type B',
                'retinal bipolar neuron type C': 'Retinal bipolar neuron type C',
                'retinal bipolar neuron type D': 'Retinal bipolar neuron type D',
                'retinal cone cell': 'Retinal cone cell',
                'retinal ganglion cell': 'Retinal ganglion cell',
                'retinal rod cell type A': 'Retinal rod cell type A',
                'retinal rod cell type B': 'Retinal rod cell type B',
                'retinal rod cell type C': 'Retinal rod cell type C',
                'unannotated': 'Unknown',
                'unspecified': 'Unknown',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "eye", "lukowski19.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                   .multiply(1/10000)
