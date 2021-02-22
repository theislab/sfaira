import anndata
import os
from typing import Union
import numpy as np

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad",
    "vieira19_Bronchi_anonymised.processed.h5ad",
]


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def __init__(
            self,
            sample_fn: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = f"human_lung_2019_10x_braga_{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_" \
                  f"10.1038/s41591-019-0468-5"

        self.download_url_data = f"https://covid19.cog.sanger.ac.uk/{self.sample_fn}"
        self.download_url_meta = None

        self.author = "Braga"
        self.doi = "10.1038/s41591-019-0468-5"
        self.healthy = True
        self.organ = "bronchus" if sample_fn == "vieira19_Bronchi_anonymised.processed.h5ad" else "lung parenchyma"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019
        self.normalization = "norm"

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

    def _load(self):
        fn = os.path.join(self.data_dir, self.sample_fn)
        adata = anndata.read(fn)
        adata.X = np.expm1(adata.X)
        self.set_unknown_class_id(ids=["1_Unicorns and artifacts"])

        return adata
