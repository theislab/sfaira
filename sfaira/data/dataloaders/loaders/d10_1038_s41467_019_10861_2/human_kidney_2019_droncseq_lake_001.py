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
        self.id = "human_kidney_2019_droncseq_lake_001_10.1038/s41467-019-10861-2"

        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121862/suppl/" \
            "GSE121862%5FUCSD%2DWU%5FSingle%5FNuclei%5FCluster%5FAnnotated%5FRaw%5FUMI%5FMatrix%2Etsv%2Egz"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121862/suppl/" \
            "GSE121862%5FUCSD%2DWU%5FSingle%5FNuclei%5FCluster%5FAnnotations%2Ecsv%2Egz"

        self.author = "Lake"
        self.doi = "10.1038/s41467-019-10861-2"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "kidney"
        self.organism = "human"
        self.protocol = "DroNc-seq"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "celltype"

        self.class_maps = {
            "0": {
                "Collecting Duct - Intercalated Cells Type A (cortex)": "Collecting Duct - Intercalated Cells Type A (cortex)",
                "Collecting Duct - Intercalated Cells Type A (medulla)": "Collecting Duct - Intercalated Cells Type A (medulla)",
                "Collecting Duct - Intercalated Cells Type B": "Collecting Duct - Intercalated Cells Type B",
                "Collecting Duct - PCs - Stressed Dissoc Subset": "Collecting Duct - PCs - Stressed Dissoc Subset",
                "Collecting Duct - Principal Cells (cortex)": "Collecting Duct - Principal Cells (cortex)",
                "Collecting Duct - Principal Cells (medulla)": "Collecting Duct - Principal Cells (medulla)",
                "Connecting Tubule": "Connecting tubule",
                "Decending Limb": "Decending Limb",
                "Distal Convoluted Tubule": "Distal Convoluted Tubule",
                "Endothelial Cells (unassigned)": "Endothelial Cells (unassigned)",
                "Endothelial Cells - AEA & DVR ": "Endothelial Cells - AEA & DVR",
                "Endothelial Cells - AVR": "Endothelial Cells - AVR",
                "Endothelial Cells - glomerular capillaries": "Endothelial Cells - glomerular capillaries",
                "Epithelial Cells (unassigned)": "Epithelial Cells (unassigned)",
                "Immune Cells - Macrophages": "Macrophage",
                "Interstitium": "Interstitium",
                "Mesangial Cells": "Mesangial Cells",
                "Podocytes": "Podocyte",
                "Proximal Tubule Epithelial Cells (S1)": "Proximal Tubule Epithelial Cells (S1)",
                "Proximal Tubule Epithelial Cells (S2)": "Proximal Tubule Epithelial Cells (S2)",
                "Proximal Tubule Epithelial Cells (S3)": "Proximal Tubule Epithelial Cells (S3)",
                "Proximal Tubule Epithelial Cells - Fibrinogen+ (S3 )": "Proximal Tubule Epithelial Cells - Fibrinogen+ (S3)",
                "Proximal Tubule Epithelial Cells - Stress/Inflam": "Proximal Tubule Epithelial Cells - Stress/Inflam",
                "Thick Ascending Limb": "Thick ascending limb of Loop of Henle",
                "Thin ascending limb": "Thin ascending limb",
                "Unknown - Novel PT CFH+ Subpopulation (S2)": "Unknown - Novel PT CFH+ Subpopulation (S2)",
                "Vascular Smooth Muscle Cells and pericytes": "Vascular Smooth Muscle Cells and pericytes",
            },
        }

    def _load(self):
        fn = [
            os.path.join(self.data_dir, "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv.gz"),
            os.path.join(self.data_dir, "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotations.csv.gz")
        ]
        adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t").T)
        annot = pd.read_csv(fn[1], index_col=0, dtype="category")
        adata.obs["celltype"] = [annot.loc[i.split("_")[0][1:]]["Annotation"] for i in adata.obs.index]

        return adata
