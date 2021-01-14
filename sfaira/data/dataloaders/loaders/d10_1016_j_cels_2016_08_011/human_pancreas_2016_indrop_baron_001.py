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
        self.id = "human_pancreas_2016_indrop_baron_001_10.1016/j.cels.2016.08.011"

        self.download = "https://covid19.cog.sanger.ac.uk/baron16.processed.h5ad"
        self.download_meta = None

        self.author = "Yanai"
        self.doi = "10.1016/j.cels.2016.08.011"
        self.healthy = True
        self.normalization = 'raw'
        self.organ = "pancreas"
        self.organism = "human"
        self.protocol = 'inDrop'
        self.state_exact = 'healthy'
        self.sub_tissue = "pancreas"
        self.year = 2016

        self.var_symbol_col = 'index'

        self.obs_key_cellontology_original = 'CellType'

        self.class_maps = {
            "0": {
                't_cell': 'T cell',
                'quiescent_stellate': 'Quiescent Stellate cell',
                'mast': 'Mast cell',
                'delta': 'Delta cell',
                'beta': 'Beta cell',
                'endothelial': 'Endothelial cell',
                'macrophage': 'Macrophage',
                'epsilon': 'Epsilon cell',
                'activated_stellate': 'Activated Stellate cell',
                'acinar': 'Acinar cell',
                'alpha': 'Alpha cell',
                'ductal': 'Ductal cell',
                'schwann': 'Schwann cell',
                'gamma': 'Gamma cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "pancreas", "baron16.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                   .multiply(1 / 10000)
