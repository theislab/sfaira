import anndata
import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = [
            "https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Mature_Full_v2.1.h5ad",
            "https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Fetal_full.h5ad"
        ]
        self.download_url_meta = None

        self.assay_sc = "10x technology"
        self.author = "Stewart"
        self.disease = "healthy"
        self.doi = "10.1126/science.aat5031"
        self.normalization = "norm"
        self.organ = "kidney"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2019

        self.gene_id_symbols_var_key = "index"
        self.gene_id_ensembl_var_key = "ID"
        self.cell_types_original_obs_key = "celltype"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "Mature_Full_v2.1.h5ad"),
        os.path.join(data_dir, "Fetal_full.h5ad")
    ]
    adult = anndata.read(fn[0])
    fetal = anndata.read(fn[1])
    adult.obs["development"] = "adult"
    fetal.obs["development"] = "fetal"
    adata = adult.concatenate(fetal)
    adata.X = np.expm1(adata.X)

    return adata
