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

        super().__init__(
            sample_id=sample_id,
            data_path=data_path,
            meta_path=meta_path,
            cache_path=cache_path,
            **kwargs
        )

        self.obs_key_sample = "derived_organ_parts_label"
        self.id = f"human_{'blood' if sample_id == 'umbilical cord blood' else 'bone'}_2018_10x_ica_" \
                  f"{str(SAMPLE_IDS.index(self.sample_id)).zfill(3)}_unknown"

        self.download_url_data = "https://data.humancellatlas.org/project-assets/project-matrices/" \
                                 "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom"
        self.download_url_meta = None

        self.author = "Regev"
        self.doi = "no_doi"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "blood" if sample_id == "umbilical cord blood" else "bone marrow"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "index"
        self.var_ensembl_col = "Accession"

        self.class_maps = {
            "0": {},
        }

    def _load_full(self):
        fn = os.path.join(self.data_dir, "cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.loom")
        self.adata = anndata.read_loom(fn)
        self.adata = self.adata[self.adata.obs["emptydrops_is_cell"] == "t"].copy()
