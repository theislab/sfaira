import anndata
import os
from typing import Union
import scipy.sparse

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "madissoon19_lung.processed.h5ad",
    "oesophagus.cellxgene.h5ad",
    "spleen.cellxgene.h5ad",
]


class Dataset(DatasetBaseGroupLoadingManyFiles):
    """
    ToDo: patient information in .obs["patient"] and sample information in .obs["sample"] (more samples than patients)
    """

    def __init__(
            self,
            sample_fn: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        organ = "lung parenchyma" if self.sample_fn == "madissoon19_lung.processed.h5ad" else \
            "esophagus" if self.sample_fn == "oesophagus.cellxgene.h5ad" else "spleen"
        self.id = f"human_{''.join(organ.split(' '))}_2019_10x_madissoon_" \
                  f"{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_10.1186/s13059-019-1906-x"

        if self.sample_fn == "madissoon19_lung.processed.h5ad":
            self.download_url_data = "https://covid19.cog.sanger.ac.uk/madissoon19_lung.processed.h5ad"
            self.var_ensembl_col = "gene.ids.HCATisStab7509734"
        elif self.sample_fn == "oesophagus.cellxgene.h5ad":
            self.download_url_data = \
                "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/oesophagus.cellxgene.h5ad"
            # Associated DCP: https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7
            self.var_ensembl_col = "gene_ids-HCATisStab7413619"
        else:
            self.download_url_data = \
                "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/spleen.cellxgene.h5ad"
            self.var_ensembl_col = "gene_ids-HCATisStab7463846"

        self.download_url_meta = None
        self.author = "Madissoon"
        self.doi = "10.1186/s13059-019-1906-x"
        self.healthy = True
        self.normalization = "raw"  # ToDo "madissoon19_lung.processed.h5ad" is close to integer but not quire (~1e-4)
        self.organ = organ
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "Celltypes"

    def _load(self):
        fn = os.path.join(self.data_dir, self.sample_fn)
        adata = anndata.read(fn)
        if self.sample_fn != "madissoon19_lung.processed.h5ad":
            adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None]))\
                .multiply(1 / 10000)
        # Cell type column called differently in madissoon19_lung.processed.h5ad:
        if self.sample_fn == "madissoon19_lung.processed.h5ad":
            adata.obs["Celltypes"] = adata.obs["CellType"]
            del adata.obs["CellType"]

        self.set_unknown_class_id(ids=["B_T_doublet", "CD34_progenitor", "Stroma"])

        return adata
