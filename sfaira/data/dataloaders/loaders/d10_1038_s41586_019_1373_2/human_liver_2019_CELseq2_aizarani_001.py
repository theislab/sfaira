import anndata
import os
from typing import Union
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_liver_2019_mCELSeq2_aizarani_001_10.1038/s41586-019-1373-2"

        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5FNormalhumanlivercellatlasdata%2Etxt%2Egz"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5Fclusterpartition%2Etxt%2Egz"

        self.author = "Gruen"
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

        self.class_maps = {
            "0": {
                "1": "NK, NKT and T cells",
                "2": "Kupffer Cell",
                "3": "NK, NKT and T cells",
                "4": "Cholangiocytes",
                "5": "NK, NKT and T cells",
                "6": "Kupffer Cell",
                "7": "Cholangiocytes",
                "8": "B Cell",
                "9": "Liver sinusoidal endothelial cells",
                "10": "Macrovascular endothelial cells",
                "11": "Hepatocyte",
                "12": "NK, NKT and T cells",
                "13": "Liver sinusoidal endothelial cells",
                "14": "Hepatocyte",
                "15": "Other endothelial cells",
                "16": "Unknown",
                "17": "Hepatocyte",
                "18": "NK, NKT and T cells",
                "19": "Unknown",
                "20": "Liver sinusoidal endothelial cells",
                "21": "Macrovascular endothelial cells",
                "22": "B Cell",
                "23": "Kupffer Cell",
                "24": "Cholangiocytes",
                "25": "Kupffer Cell",
                "26": "Other endothelial cells",
                "27": "Unknown",
                "28": "NK, NKT and T cells",
                "29": "Macrovascular endothelial cells",
                "30": "Hepatocyte",
                "31": "Kupffer Cell",
                "32": "Liver sinusoidal endothelial cells",
                "33": "Hepatic stellate cells",
                "34": "B Cell",
                "35": "Other endothelial cells",
                "36": "Unknown",
                "37": "Unknown",
                "38": "B Cell",
                "39": "Cholangiocytes"
            },
        }

    def _load(self):
        fn = [
            os.path.join(self.doi_path, "GSE124395_Normalhumanlivercellatlasdata.txt.gz"),
            os.path.join(self.doi_path, "GSE124395_clusterpartition.txt.gz")
        ]
        self.adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t").T)
        celltype_df = pd.read_csv(fn[1], sep=" ")
        self.adata = self.adata[[i in celltype_df.index for i in self.adata.obs.index]].copy()
        self.adata.obs["CellType"] = [str(celltype_df.loc[i]["sct@cpart"]) for i in self.adata.obs.index]
