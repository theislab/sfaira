import anndata
import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/voigt19.processed.h5ad"
        self.download_url_meta = None

        self.author = "Voigt"
        self.doi = "10.1073/pnas.1914143116"
        self.healthy = True
        self.normalization = "norm"
        self.organ = "retina"
        self.organism = "human"
        self.assay_sc = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"
        self.cellontology_original_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "voigt19.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)

    return adata
