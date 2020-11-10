import os
from typing import Union
from .external import DatasetBase
import anndata
import numpy as np
import scipy.sparse


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
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_pancreas_2016_indrop_baron_001_10.1016/j.cels.2016.08.011"
        self.download_website = "https://covid19.cog.sanger.ac.uk/baron16.processed.h5ad"
        self.organ = "pancreas"
        self.sub_tissue = "pancreas"
        self.has_celltypes = True

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
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "pancreas", "baron16.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns["lab"] = "Yanai"
        self.adata.uns["year"] = 2016
        self.adata.uns["doi"] = "10.1016/j.cels.2016.08.011"
        self.adata.uns["protocol"] = 'inDrop'
        self.adata.uns["organ"] = self.organ
        self.adata.uns["subtissue"] = self.sub_tissue
        self.adata.uns["animal"] = "human"
        self.adata.uns["id"] = self.id
        self.adata.uns["wget_download"] = self.download_website
        self.adata.uns["has_celltypes"] = self.has_celltypes
        self.adata.uns["counts"] = 'raw'
        self.adata.obs["healthy"] = True
        self.adata.obs["state_exact"] = 'healthy'

        self.adata.obs["cell_ontology_class"] = self.adata.obs['CellType']
        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index='ensembl')
