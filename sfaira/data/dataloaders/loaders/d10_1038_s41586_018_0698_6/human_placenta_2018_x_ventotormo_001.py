import os
from typing import Union
import pandas as pd
import anndata

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "E-MTAB-6678.processed",
    "E-MTAB-6701.processed",
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
        self.download_url_data = f"https://www.ebi.ac.uk/arrayexpress/files/{self.sample_fn.split('.')[0]}/" \
                                 f"{self.sample_fn}.1.zip"
        self.download_url_meta = f"https://www.ebi.ac.uk/arrayexpress/files/{self.sample_fn.split('.')[0]}/" \
                                 f"{self.sample_fn}.2.zip"

        self.author = "Ventotormo"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "placenta"
        self.organism = "human"
        self.doi = "10.1038/s41586-018-0698-6"
        self.protocol = "10X sequencing" if self.sample_fn == "E-MTAB-6678.processed" else "Smart-seq2"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "names"
        self.var_ensembl_col = "ensembl"
        self.obs_key_cellontology_original = "annotation"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = [
        os.path.join(data_dir, f"{sample_fn}.1.zip"),
        os.path.join(data_dir, f"{sample_fn}.2.zip"),
    ]
    adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t", index_col="Gene").T)
    df = pd.read_csv(fn[1], sep="\t")
    for i in df.columns:
        adata.obs[i] = [df.loc[j][i] for j in adata.obs.index]

    adata.var["ensembl"] = [i.split("_")[1] for i in adata.var.index]
    adata.var["names"] = [i.split("_")[0] for i in adata.var.index]
    adata.var = adata.var.reset_index().reset_index().drop("index", axis=1)
    adata = adata[:, ~adata.var.index.isin(
        ["", "-1", "-10", "-11", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "A.2", "A.3"])].copy()

    return adata
