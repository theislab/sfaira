import anndata
import os
import numpy as np

from sfaira.data import DatasetBase

SAMPLE_FNS = [
    "vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad",
    "vieira19_Bronchi_anonymised.processed.h5ad",
]


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = f"https://covid19.cog.sanger.ac.uk/{self.sample_fn}"
        self.download_url_meta = None

        self.author = "Braga"
        self.doi = "10.1038/s41591-019-0468-5"
        self.healthy = True
        self.organ = "bronchus" if self.sample_fn == "vieira19_Bronchi_anonymised.processed.h5ad" else "lung parenchyma"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019
        self.normalization = "norm"

        self.var_symbol_col = "index"
        self.cellontology_original_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)

    return adata
