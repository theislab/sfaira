import anndata
import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    """
    TODO transform field development to controlled field age
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = [
            "https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Mature_Full_v2.1.h5ad",
            "https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Fetal_full.h5ad"
        ]
        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Stewart"
        self.disease = "healthy"
        self.doi_journal = "10.1126/science.aat5031"
        self.layer_processed = "X"
        self.organ = "kidney"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2019

        self.feature_id_var_key = "ID"
        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"
        self.cell_type_obs_key = "celltype"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "Mature_Full_v2.1.h5ad"),
        os.path.join(data_dir, "Fetal_full.h5ad")
    ]
    adult = anndata.read(fn[0])
    fetal = anndata.read(fn[1])
    # TODO this is is not a controlled field
    adult.obs["development"] = "adult"
    fetal.obs["development"] = "fetal"
    adata = adult.concatenate(fetal)
    adata.X = np.expm1(adata.X)

    return adata
