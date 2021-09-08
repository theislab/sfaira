import anndata
import os
import scipy.sparse

from sfaira.data import DatasetBase

SAMPLE_FNS = [
    "madissoon19_lung.processed.h5ad",
    "oesophagus.cellxgene.h5ad",
    "spleen.cellxgene.h5ad",
]


class Dataset(DatasetBase):
    """
    ToDo: patient information in .obs["patient"] and sample information in .obs["sample"] (more samples than patients)
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.sample_fn == "madissoon19_lung.processed.h5ad":
            self.download_url_data = "https://covid19.cog.sanger.ac.uk/madissoon19_lung.processed.h5ad"
            self.gene_id_ensembl_var_key = "gene.ids.HCATisStab7509734"
        elif self.sample_fn == "oesophagus.cellxgene.h5ad":
            self.download_url_data = \
                "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/oesophagus.cellxgene.h5ad"
            # Associated DCP: https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7
            self.gene_id_ensembl_var_key = "gene_ids-HCATisStab7413619"
        else:
            self.download_url_data = \
                "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/spleen.cellxgene.h5ad"
            self.gene_id_ensembl_var_key = "gene_ids-HCATisStab7463846"

        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Madissoon"
        self.disease = "healthy"
        self.doi_journal = "10.1186/s13059-019-1906-x"
        self.normalization = "raw"  # ToDo "madissoon19_lung.processed.h5ad" is close to integer but not quire (~1e-4)
        self.organ = "lung parenchyma" if self.sample_fn == "madissoon19_lung.processed.h5ad" else \
            "esophagus" if self.sample_fn == "oesophagus.cellxgene.h5ad" else "spleen"
        self.organism = "human"
        self.year = 2019
        self.sample_source = "primary_tissue"

        self.gene_id_symbols_var_key = "index"
        self.cell_type_obs_key = "Celltypes"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    if sample_fn != "madissoon19_lung.processed.h5ad":
        adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None]))\
            .multiply(1 / 10000)
    # Cell type column called differently in madissoon19_lung.processed.h5ad:
    if sample_fn == "madissoon19_lung.processed.h5ad":
        adata.obs["Celltypes"] = adata.obs["CellType"]
        del adata.obs["CellType"]

    return adata
