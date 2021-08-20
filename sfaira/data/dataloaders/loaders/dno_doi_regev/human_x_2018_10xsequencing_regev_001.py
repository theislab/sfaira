import anndata
import os

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    """
    TODO data link is outdated. Maybe update to this
    https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79/project-matrices.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://data.humancellatlas.org/project-assets/project-matrices/" \
                                 "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Regev"
        self.disease = "healthy"
        self.doi_journal = "no_doi_regev"
        self.normalization = "raw"
        self.organ_obs_key = "derived_organ_parts_label"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.year = 2018

        self.gene_id_symbols_var_key = "index"
        self.gene_id_ensembl_var_key = "Accession"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom")
    adata = anndata.read_loom(fn)
    adata = adata[adata.obs["emptydrops_is_cell"] == "t"].copy()

    return adata
