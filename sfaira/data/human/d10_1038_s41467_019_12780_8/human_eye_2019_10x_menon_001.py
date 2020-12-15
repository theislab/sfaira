import anndata
import os
from typing import Union
from .external import DatasetBase


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
        super().__init__(path=path, meta_path=meta_path, **kwargs)
        self.organism = "human"
        self.id = "human_eye_2019_10x_menon_001_10.1038/s41467-019-12780-8"
        self.download = "https://covid19.cog.sanger.ac.uk/menon19.processed.h5ad"
        self.download_meta = None
        self.organ = "eye"
        self.sub_tissue = "retina"
        self.author = 'Hafler'
        self.year = 2019
        self.doi = '10.1038/s41467-019-12780-8'
        self.protocol = '10x'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'CellType'

        self.class_maps = {
            "0": {
                'ACs': 'Amacrine cell',
                'BPs': 'BPs',
                'Cones': 'Retinal cone cell',
                'Endo': 'Endothelial cell',
                'HCs': 'Horizontal cells',
                'Macroglia': 'Macroglia',
                'Microglia': 'Microglia',
                'RGCs': 'Retinal ganglion cell',
                'Rods': 'Rods',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "eye", "menon19.processed.h5ad")
        self.adata = anndata.read(fn)
