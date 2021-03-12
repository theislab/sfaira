import anndata
import os

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://data.humancellatlas.org/project-assets/project-matrices/" \
                                 "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom"
        self.download_url_meta = None

        self.author = "Regev"
        self.doi = "no_doi_regev"
        self.healthy = True
        self.normalization = "raw"
        self.organ_obs_key = "derived_organ_parts_label"
        self.organism = "human"
        self.assay_sc = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2018
        self.sample_source = "primary_tissue"

        self.var_symbol_col = "index"
        self.var_ensembl_col = "Accession"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom")
    adata = anndata.read_loom(fn)
    adata = adata[adata.obs["emptydrops_is_cell"] == "t"].copy()

    return adata
