import anndata
import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/park20.processed.h5ad"
        self.download_url_meta = None

        self.author = "Park"
        self.doi = "10.1126/science.aay3224"
        self.healthy = True
        self.normalization = "norm"
        self.organ = "thymus"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"
        self.cellontology_original_obs_key = "Anno_level_fig1"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "park20.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)

    return adata
