import anndata
import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/voigt19.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v3"
        self.author = "Voigt"
        self.disease = "healthy"
        self.doi_journal = "10.1073/pnas.1914143116"
        self.layer_processed = "X"
        self.organ = "retina"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.year = 2019

        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"
        self.cell_type_obs_key = "CellType"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "voigt19.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)

    return adata
