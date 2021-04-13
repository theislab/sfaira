import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    ToDo: revisit gamma cell missing in CO
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/baron16.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc = "inDrop"
        self.author = "Baron"
        self.disease = "healthy"
        self.doi = "10.1016/j.cels.2016.08.011"
        self.normalization = "raw"
        self.organ = "pancreas"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2016

        self.gene_id_symbols_var_key = "index"
        self.cell_types_original_obs_key = "CellType"

        self.set_dataset_id(idx=1)

        self.class_maps = {
            "0": {
                "t_cell": "T cell",
                "quiescent_stellate": "Quiescent Stellate cell",
                "mast": "Mast cell",
                "delta": "Delta cell",
                "beta": "Beta cell",
                "endothelial": "Endothelial cell",
                "macrophage": "Macrophage",
                "epsilon": "Epsilon cell",
                "activated_stellate": "Activated Stellate cell",
                "acinar": "Acinar cell",
                "alpha": "Alpha cell",
                "ductal": "Ductal cell",
                "schwann": "Schwann cell",
                "gamma": "Gamma cell",
            },
        }


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "baron16.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)

    return adata
