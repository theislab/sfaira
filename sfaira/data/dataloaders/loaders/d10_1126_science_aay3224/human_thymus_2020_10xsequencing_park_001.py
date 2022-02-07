import anndata
import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    """
    TODO add meta data

    .obs columns Age contains entries ['3m', '6m', '7w', '8w', '9w', '10m', '10w', '11w', '12w', '13w', '13y',
       '14w', '15m', '16w', '17w', '24y', '30m', '35y']
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/park20.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc = "10x technology"
        self.author = "Park"
        self.disease = "healthy"
        self.doi_journal = "10.1126/science.aay3224"
        self.individual_obs_key = "donor"
        self.layer_processed = "X"
        self.organ = "thymus"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.sex_obs_key = "Gender"
        self.tech_sample_obs_key = "Sample"
        self.year = 2020

        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"
        self.cell_type_obs_key = "Anno_level_fig1"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "park20.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.obs["Gender"] = [{"Male": "male", "Female": "female"}[x] for x in adata.obs["Gender"].values]

    return adata
