import anndata
import os
from typing import Union
import numpy as np

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
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/voigt19.processed.h5ad"
        self.download_url_meta = None

        self.author = "Voigt"
        self.doi = "10.1073/pnas.1914143116"
        self.healthy = True
        self.normalization = "norm"
        self.organ = "retina"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "B-cell": "B-cell",
                "Endothelial": "Endothelial cell",
                "Fibroblast": "Fibroblast",
                "Macrophage": "Macrophage",
                "Mast-cell": "Mast-cell",
                "Melanocyte": "Melanocyte",
                "Pericyte": "Pericyte",
                "RPE": "Retinal pigment epithelium",
                "Schwann1": "Schwann1",
                "Schwann2": "Schwann2",
                "T/NK-cell": "T/NK-cell",
            },
        }

    def _load(self):
        fn = os.path.join(self.data_dir, "voigt19.processed.h5ad")
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)

        return adata
