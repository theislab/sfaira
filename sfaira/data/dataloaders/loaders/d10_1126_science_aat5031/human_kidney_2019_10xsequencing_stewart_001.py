import anndata
import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    """
    TODO: annotate developmental cell types in set_unknown_class_id
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = [
            "https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Mature_Full_v2.1.h5ad",
            "https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Fetal_full.h5ad"
        ]
        self.download_url_meta = None

        self.author = "Stewart"
        self.doi = "10.1126/science.aat5031"
        self.healthy = True
        self.normalization = "norm"
        self.organ = "kidney"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"
        self.var_ensembl_col = "ID"
        self.obs_key_cellontology_original = "celltype"

        self.set_dataset_id(idx=1)

        self.set_unknown_class_id(ids=[
            "CNT/PC - proximal UB", "Distal S shaped body", "Medial S shaped body", "Proliferating stroma progenitor",
            "Proximal S shaped body", "Stroma progenitor", "Proximal UB",
        ])


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
