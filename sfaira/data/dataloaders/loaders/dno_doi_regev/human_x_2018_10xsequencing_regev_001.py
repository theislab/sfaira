import anndata
import os
from typing import Union

from sfaira.data import DatasetBaseGroupLoadingOneFile

SAMPLE_IDS = [
    "umbilical cord blood",
    "bone marrow"
]


class Dataset(DatasetBaseGroupLoadingOneFile):

    def __init__(
            self,
            sample_id: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_id=sample_id, data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.obs_key_sample = "derived_organ_parts_label"

        self.download_url_data = "https://data.humancellatlas.org/project-assets/project-matrices/" \
                                 "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom"
        self.download_url_meta = None

        self.author = "Regev"
        self.doi = "no_doi_regev"
        self.healthy = True
        self.normalization = "raw"
        self.organ = sample_id
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "index"
        self.var_ensembl_col = "Accession"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom")
    adata = anndata.read_loom(fn)
    adata = adata[adata.obs["emptydrops_is_cell"] == "t"].copy()

    return adata
