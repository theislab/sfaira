import anndata
import os

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/menon19.processed.h5ad"
        self.download_url_meta = None

        self.author = "Menon"
        self.doi = "10.1038/s41467-019-12780-8"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "retina"
        self.organism = "human"
        self.assay_sc = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"
        self.cellontology_original_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "menon19.processed.h5ad")
    adata = anndata.read(fn)

    return adata
