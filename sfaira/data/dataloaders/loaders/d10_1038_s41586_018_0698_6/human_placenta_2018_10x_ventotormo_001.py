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
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        protocol = "10x" if self.sample_fn == "E-MTAB-6678.processed" else "smartseq2"
        self.id = f"human_placenta_2018_{protocol}_ventotormo_{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_" \
                  f"10.1038/s41586-018-0698-6"

        self.download_url_data = f"https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6701/{self.sample_fn}.1.zip"
        self.download_url_meta = f"https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6701/{self.sample_fn}.2.zip"

        self.author = "Teichmann"
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

        self.class_maps = {
            "0": {
                "DC1": "Dendritic Cells 1",
                "DC2": "Dendritic Cells 2",
                "EVT": "Extravillous Trophoblasts",
                "Endo (f)": "Endothelial Cells f",
                "Endo (m)": "Endothelial Cells m",
                "Endo L": "Endothelial Cells L",
                "Epi1": "Epithelial Glandular Cells 1",
                "Epi2": "Epithelial Glandular Cells 2",
                "Granulocytes": "Granulocytes",
                "HB": "Hofbauer Cells",
                "ILC3": "ILC3",
                "MO": "Monocyte",
                "NK CD16+": "NK Cells CD16+",
                "NK CD16-": "NK Cells CD16-",
                "Plasma": "B cell (Plasmocyte)",
                "SCT": "Syncytiotrophoblasts",
                "Tcells": "T cell",
                "VCT": "Villous Cytotrophoblasts",
                "dM1": "Decidual Macrophages 1",
                "dM2": "Decidual Macrophages 2",
                "dM3": "Decidual Macrophages 3",
                "dNK p": "Decidual NK Cells p",
                "dNK1": "Decidual NK Cells 1",
                "dNK2": "Decidual NK Cells 2",
                "dNK3": "Decidual NK Cells 3",
                "dP1": "Perivascular Cells 1",
                "dP2": "Perivascular Cells 2",
                "dS1": "Decidual Stromal Cells 1",
                "dS2": "Decidual Stromal Cells 2",
                "dS3": "Decidual Stromal Cells 3",
                "fFB1": "Fibroblasts 1",
                "fFB2": "Fibroblasts 2",
            },
        }

    def _load(self):
        fn = [
            os.path.join(self.doi_path, f"{self.sample_fn}.1.zip"),
            os.path.join(self.doi_path, f"{self.sample_fn}.2.zip"),
        ]
        self.adata = anndata.AnnData(pd.read_csv(fn[0], sep="\t", index_col="Gene").T)
        df = pd.read_csv(fn[1], sep="\t")
        for i in df.columns:
            self.adata.obs[i] = [df.loc[j][i] for j in self.adata.obs.index]

        self.adata.var["ensembl"] = [i.split("_")[1] for i in self.adata.var.index]
        self.adata.var["names"] = [i.split("_")[0] for i in self.adata.var.index]
        self.adata.var = self.adata.var.reset_index().reset_index().drop("index", axis=1)
        self.adata = self.adata[:, ~self.adata.var.index.isin(
            ["", "-1", "-10", "-11", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "A.2", "A.3"])].copy()
