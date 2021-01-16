import anndata
import os
from typing import Union
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    The input files for this dataloader (GSE115469.csv.gz and GSE115469_labels.txt) were kindly provided to us by the
    authors of the publication. Please contact them directly to obtain the required
    files.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_liver_2018_10x_macparland_001_10.1038/s41467-018-06318-7"

        self.download = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469"
        self.download_meta = "private"

        self.author = "McGilvray"
        self.doi = "10.1038/s41467-018-06318-7"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "liver"  # ToDo: "caudate lobe"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "celltype"

        self.class_maps = {
            "0": {
                "1": "Hepatocyte 1",
                "2": "Alpha beta T cells",
                "3": "Hepatocyte 2",
                "4": "Inflammatory macrophages",
                "5": "Hepatocyte 3",
                "6": "Hepatocyte 4",
                "7": "Plasma cells",
                "8": "NK cell",
                "9": "Gamma delta T cells 1",
                "10": "Non inflammatory macrophages",
                "11": "Periportal LSECs",
                "12": "Central venous LSECs",
                "13": "Endothelial cell",
                "14": "Hepatocyte 5",
                "15": "Hepatocyte 6",
                "16": "Mature B cells",
                "17": "Cholangiocytes",
                "18": "Gamma delta T cells 2",
                "19": "Erythroid cells",
                "20": "Hepatic stellate cells"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = [
                os.path.join(self.path, "human", "liver", "GSE115469.csv.gz"),
                os.path.join(self.path, "human", "liver", "GSE115469_labels.txt")
            ]
        self.adata = anndata.read_csv(fn[0]).T
        celltype_df = pd.read_csv(fn[1], sep="\t").set_index("CellName")
        self.adata.obs["celltype"] = [str(celltype_df.loc[i]["Cluster#"]) for i in self.adata.obs.index]
