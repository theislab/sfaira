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
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5FNormalhumanlivercellatlasdata%2Etxt%2Egz"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5Fclusterpartition%2Etxt%2Egz"

        self.author = "Aizarani"
        self.doi = "10.1038/s41586-019-1373-2"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "liver"
        self.organism = "human"
        self.protocol = "CEL-seq2"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"
        self.obs_key_cellontology_original = "CellType"

        self.set_dataset_id(idx=1)

    def _load(self):
        fn = [
            os.path.join(self.data_dir, "GSE124395_Normalhumanlivercellatlasdata.txt.gz"),
            os.path.join(self.data_dir, "GSE124395_clusterpartition.txt.gz")
        ]
        adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t").T)
        celltype_df = pd.read_csv(fn[1], sep=" ")
        adata = adata[[i in celltype_df.index for i in adata.obs.index]].copy()
        adata.obs["CellType"] = [str(celltype_df.loc[i]["sct@cpart"]) for i in adata.obs.index]

        self.set_unknown_class_id(ids=["16", "19", "27", "36", "37"])

        return adata
