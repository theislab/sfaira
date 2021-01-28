import anndata
import os
from typing import Union
import tarfile
import pandas as pd
import scipy.sparse

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "GSM3589406_PP001swap.filtered.matrix.txt.gz",
    "GSM3589407_PP002swap.filtered.matrix.txt.gz",
    "GSM3589408_PP003swap.filtered.matrix.txt.gz",
    "GSM3589409_PP004swap.filtered.matrix.txt.gz",
    "GSM3589410_PP005swap.filtered.matrix.txt.gz",
    "GSM3589411_PP006swap.filtered.matrix.txt.gz",
    "GSM3589412_PP009swap.filtered.matrix.txt.gz",
    "GSM3589413_PP010swap.filtered.matrix.txt.gz",
    "GSM3589414_PP011swap.filtered.matrix.txt.gz",
    "GSM3589415_PP012swap.filtered.matrix.txt.gz",
    "GSM3589416_PP013swap.filtered.matrix.txt.gz",
    "GSM3589417_PP014swap.filtered.matrix.txt.gz",
    "GSM3589418_PP017swap.filtered.matrix.txt.gz",
    "GSM3589419_PP018swap.filtered.matrix.txt.gz",
    "GSM3589420_PP019swap.filtered.matrix.txt.gz",
    "GSM3589421_PP020swap.filtered.matrix.txt.gz",
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
        self.id = f"human_mixed_2019_10x_szabo_{str(SAMPLE_FNS.index(sample_fn)).zfill(3)}_10.1038/s41467-019-12464-3"

        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126030/suppl/GSE126030_RAW.tar"
        self.download_url_meta = [
            "private,donor1.annotation.txt",
            "private,donor2.annotation.txt"
        ]

        self.author = "Sims"
        self.doi = "10.1038/s41467-019-12464-3"
        self.healthy = True
        self.normalization = "raw"
        self.organism = "human"
        self.protocol = "10x"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "Gene"
        self.var_ensembl_col = "Accession"

        self.obs_key_cellontology_original = "cell_ontology_class"
        self.obs_key_organ = "organ"

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        fn_tar = os.path.join(self.path, "human", "mixed", "GSE126030_RAW.tar"),
        fn_annot = [
            os.path.join(self.path, "human", "mixed", "donor1.annotation.txt"),
            os.path.join(self.path, "human", "mixed", "donor2.annotation.txt")
        ]
        with tarfile.open(fn_tar) as tar:
            df = pd.read_csv(tar.extractfile(self.sample_fn), compression="gzip", sep="\t")
            df.index = [i.split(".")[0] for i in df["Accession"]]
            var = pd.concat([df.pop(x) for x in ["Gene", "Accession"]], 1)
            if df.columns[-1].startswith("Un"):
                df.drop(df.columns[len(df.columns) - 1], axis=1, inplace=True)
            self.adata = anndata.AnnData(df.T)
            self.adata.var = var
            if "PP001" in fn or "PP002" in fn:
                self.adata.obs["donor"] = "Donor1"
                self.adata.obs["organ"] = "Lung"
            elif "PP003" in fn or "PP004" in fn:
                self.adata.obs["donor"] = "Donor1"
                self.adata.obs["organ"] = "Bone Marrow"
            elif "PP005" in fn or "PP006" in fn:
                self.adata.obs["donor"] = "Donor1"
                self.adata.obs["organ"] = "Lymph Node"
            elif "PP009" in fn or "PP010" in fn:
                self.adata.obs["donor"] = "Donor2"
                self.adata.obs["organ"] = "Lung"
            elif "PP011" in fn or "PP012" in fn:
                self.adata.obs["donor"] = "Donor2"
                self.adata.obs["organ"] = "Bone Marrow"
            elif "PP013" in fn or "PP014" in fn:
                self.adata.obs["donor"] = "Donor2"
                self.adata.obs["organ"] = "Lymph Node"
            self.adata.obs.index = fn.split("_")[1].split("s")[0] + "nskept." + self.adata.obs.index
        self.adata.obs["cell_ontology_class"] = "Unknown"
        df1 = pd.read_csv(fn_annot[0], sep="\t", index_col=0, header=None)
        df2 = pd.read_csv(fn_annot[1], sep="\t", index_col=0, header=None)
        for i in df1.index:
            self.adata.obs["cell_ontology_class"].loc[i] = df1.loc[i][1]
        for i in df2.index:
            self.adata.obs["cell_ontology_class"].loc[i] = df2.loc[i][1]
        self.adata.X = scipy.sparse.csc_matrix(self.adata.X)
