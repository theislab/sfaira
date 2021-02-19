import anndata
import os
from typing import Union
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_liver_2018_10xsequencing_macparland_001_10.1038/s41467-018-06318-7"

        self.download_url_data = "private,GSE115469.csv.gz"
        self.download_url_meta = "private,GSE115469_labels.txt"

        self.author = "MacParland"
        self.doi = "10.1038/s41467-018-06318-7"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "caudate lobe of liver"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "celltype"

    def _load(self):
        fn = [
            os.path.join(self.data_dir, "GSE115469.csv.gz"),
            os.path.join(self.data_dir, "GSE115469_labels.txt")
        ]
        adata = anndata.read_csv(fn[0]).T
        celltype_df = pd.read_csv(fn[1], sep="\t").set_index("CellName")
        adata.obs["celltype"] = [str(celltype_df.loc[i]["Cluster#"]) for i in adata.obs.index]

        return adata
