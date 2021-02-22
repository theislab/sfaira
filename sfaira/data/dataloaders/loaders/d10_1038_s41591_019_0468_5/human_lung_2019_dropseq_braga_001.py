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
        self.id = "human_lung_2019_dropseq_braga_003_10.1038/s41591-019-0468-5"

        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fraw%5Fcounts%2Ecsv%2Egz"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fbarcodes%5Fcell%5Ftypes%2Etxt%2Egz"

        self.author = "Braga"
        self.doi = "10.1038/s41591-019-0468-5"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.protocol = "Drop-seq"
        self.state_exact = "uninvolved areas of tumour resection material"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "celltype"

    def _load(self):
        fn = [
            os.path.join(self.data_dir, "GSE130148_raw_counts.csv.gz"),
            os.path.join(self.data_dir, "GSE130148_barcodes_cell_types.txt.gz"),
        ]
        adata = anndata.read_csv(fn[0]).T
        adata.obs = pd.read_csv(fn[1], sep="\t", index_col=0)
        self.set_unknown_class_id(ids=["1_Unicorns and artifacts"])

        return adata
