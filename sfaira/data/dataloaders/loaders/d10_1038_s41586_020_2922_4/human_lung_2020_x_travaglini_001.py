import anndata
import os
import scipy.sparse
import numpy as np

from sfaira.data import DatasetBase

SAMPLE_FNS = [
    "droplet_normal_lung_blood_scanpy.20200205.RC4.h5ad",
    "facs_normal_lung_blood_scanpy.20200205.RC4.h5ad"
]


class Dataset(DatasetBase):
    """
    ToDo split by sample / patient in obs columns:
      bio replicates droplet file "orig.ident"+"sample"+"magnetic.selection",
      bio replicates facs file "patient"+"sample"
      tech replicates droplet file "channel",
      tech replicates facs file "plate.barcode"
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        synapse_id = {
            "droplet_normal_lung_blood_scanpy.20200205.RC4.h5ad": "syn21625095",
            "facs_normal_lung_blood_scanpy.20200205.RC4.h5ad": "syn21625142"
        }

        self.download_url_data = f"{synapse_id[self.sample_fn]},{self.sample_fn}"
        self.download_url_meta = None

        self.author = "Travaglini"
        self.doi = "10.1038/s41586-020-2922-4"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.protocol = "10X sequencing" if self.sample_fn.split("_")[0] == "droplet" else "Smart-seq2"
        self.state_exact = "healthy"
        self.year = 2020

        self.obs_key_cellontology_original = "free_annotation"
        self.var_symbol_col = "index"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    if sample_fn.split("_")[0] == "droplet":
        norm_const = 10000
        sf_key = "nUMI"
    else:
        norm_const = 1000000
        sf_key = "nReads"
    adata = anndata.read(fn)
    adata.X = scipy.sparse.csc_matrix(adata.X)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs[sf_key].values[:, None])).multiply(1 / norm_const)

    return adata
