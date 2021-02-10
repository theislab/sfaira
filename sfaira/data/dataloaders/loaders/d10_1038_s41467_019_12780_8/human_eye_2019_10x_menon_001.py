import anndata
import os
from typing import Union

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_eye_2019_10x_menon_001_10.1038/s41467-019-12780-8"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/menon19.processed.h5ad"
        self.download_url_meta = None

        self.author = "Hafler"
        self.doi = "10.1038/s41467-019-12780-8"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "eye"  # ToDo: "retina"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "ACs": "Amacrine cell",
                "BPs": "BPs",
                "Cones": "Retinal cone cell",
                "Endo": "Endothelial cell",
                "HCs": "Horizontal cells",
                "Macroglia": "Macroglia",
                "Microglia": "Microglia",
                "RGCs": "Retinal ganglion cell",
                "Rods": "Rods",
            },
        }

    def _load(self):
        fn = os.path.join(self.full_path, "menon19.processed.h5ad")
        self.adata = anndata.read(fn)
