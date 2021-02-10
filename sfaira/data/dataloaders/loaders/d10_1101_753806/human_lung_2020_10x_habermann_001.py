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
        self.id = "human_lung_2020_10x_habermann_001_10.1101/753806"

        self.download_url_data = [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fmatrix%2Emtx%2Egz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fgenes%2Etsv%2Egz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fbarcodes%2Etsv%2Egz"
        ]
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5FIPF%5Fmetadata%2Ecsv%2Egz"

        self.author = "Kropski"
        self.doi = "10.1101/753806"
        self.normalization = "raw"
        self.organ = "lung"  # ToDo: "parenchyma"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.year = 2020

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "celltype"
        self.obs_key_state_exact = "Diagnosis"
        self.obs_key_healthy = "Status"
        self.healthy_state_healthy = "Control"

        self.class_maps = {
            "0": {
                "Proliferating Macrophages": "Macrophages",
                "Myofibroblasts": "Myofibroblasts",
                "Proliferating Epithelial Cells": "Proliferating Epithelial Cells",
                "Mesothelial Cells": "Mesothelium",
                "cDCs": "Dendritic cells",
                "Mast Cells": "Mast cells",
                "Ciliated": "Multiciliated lineage",
                "T Cells": "T cell lineage",
                "pDCs": "Dendritic cells",
                "Smooth Muscle Cells": "2_Smooth Muscle",
                "Transitional AT2": "AT2",
                "AT2": "AT2",
                "B Cells": "B cell lineage",
                "NK Cells": "Innate lymphoid cells",
                "Monocytes": "Monocytes",
                "Basal": "Basal",
                "Plasma Cells": "B cell lineage",
                "Differentiating Ciliated": "Multiciliated lineage",
                "Macrophages": "Macrophages",
                "MUC5B+": "Secretory",
                "SCGB3A2+": "Secretory",
                "Fibroblasts": "Fibroblasts",
                "Lymphatic Endothelial Cells": "Lymphatic EC",
                "Endothelial Cells": "2_Blood vessels",
                "SCGB3A2+ SCGB1A1+": "Secretory",
                "PLIN2+ Fibroblasts": "Fibroblasts",
                "KRT5-/KRT17+": "KRT5-/KRT17+",
                "MUC5AC+ High": "Secretory",
                "Proliferating T Cells": "T cell lineage",
                "AT1": "AT1",
                "HAS1 High Fibroblasts": "Fibroblasts"
            },
        }

    def _load(self):
        fn = [
            os.path.join(self.full_path, "GSE135893_matrix.mtx.gz"),
            os.path.join(self.full_path, "GSE135893_genes.tsv.gz"),
            os.path.join(self.full_path, "GSE135893_barcodes.tsv.gz"),
            os.path.join(self.full_path, "GSE135893_IPF_metadata.csv.gz"),
        ]
        self.adata = anndata.read_mtx(fn[0]).T
        self.adata.var = pd.read_csv(fn[1], index_col=0, header=None, names=["ids"])
        self.adata.obs = pd.read_csv(fn[2], index_col=0, header=None, names=["barcodes"])
        obs = pd.read_csv(fn[3], index_col=0)
        self.adata = self.adata[obs.index.tolist(), :].copy()
        self.adata.obs = obs

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
