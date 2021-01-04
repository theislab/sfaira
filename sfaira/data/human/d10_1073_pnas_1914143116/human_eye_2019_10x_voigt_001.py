import anndata
import os
from typing import Union
from .external import DatasetBase
import numpy as np


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
        self.organism = "human"
        self.id = "human_eye_2019_10x_voigt_001_10.1073/pnas.1914143116"
        self.download = "https://covid19.cog.sanger.ac.uk/voigt19.processed.h5ad"
        self.download_meta = None
        self.organ = "eye"
        self.sub_tissue = "retina"
        self.author = 'Mullins'
        self.year = 2019
        self.doi = '10.1073/pnas.1914143116'
        self.protocol = '10x'
        self.normalization = 'norm'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'CellType'

        self.class_maps = {
            "0": {
                'B-cell': 'B-cell',
                'Endothelial': 'Endothelial cell',
                'Fibroblast': 'Fibroblast',
                'Macrophage': 'Macrophage',
                'Mast-cell': 'Mast-cell',
                'Melanocyte': 'Melanocyte',
                'Pericyte': 'Pericyte',
                'RPE': 'Retinal pigment epithelium',
                'Schwann1': 'Schwann1',
                'Schwann2': 'Schwann2',
                'T/NK-cell': 'T/NK-cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "eye", "voigt19.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
