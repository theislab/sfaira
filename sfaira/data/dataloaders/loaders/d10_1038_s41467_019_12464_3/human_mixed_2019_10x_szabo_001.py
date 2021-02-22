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
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, sample_fns=SAMPLE_FNS, data_path=data_path, meta_path=meta_path,
                         cache_path=cache_path, **kwargs)
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126030/suppl/GSE126030_RAW.tar"
        self.download_url_meta = [
            "private,donor1.annotation.txt",
            "private,donor2.annotation.txt"
        ]

        self.author = "Szabo"
        self.doi = "10.1038/s41467-019-12464-3"
        self.healthy = True
        self.normalization = "raw"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_symbol_col = "Gene"
        self.var_ensembl_col = "Accession"

        self.obs_key_cellontology_original = "cell_ontology_class"
        self.obs_key_organ = "organ"

        self.class_maps = {
            "0": {},
        }

    def _load(self):
        fn = [
            os.path.join(self.data_dir, "GSE126030_RAW.tar"),
            os.path.join(self.data_dir, "donor1.annotation.txt"),
            os.path.join(self.data_dir, "donor2.annotation.txt")
        ]
        with tarfile.open(fn[0]) as tar:
            df = pd.read_csv(tar.extractfile(self.sample_fn), compression="gzip", sep="\t")
            df.index = [i.split(".")[0] for i in df["Accession"]]
            var = pd.concat([df.pop(x) for x in ["Gene", "Accession"]], 1)
            if df.columns[-1].startswith("Un"):
                df.drop(df.columns[len(df.columns) - 1], axis=1, inplace=True)
            adata = anndata.AnnData(df.T)
            adata.var = var
            if "PP001" in self.sample_fn or "PP002" in self.sample_fn:
                adata.obs["donor"] = "Donor1"
                adata.obs["organ"] = "lung"
            elif "PP003" in self.sample_fn or "PP004" in self.sample_fn:
                adata.obs["donor"] = "Donor1"
                adata.obs["organ"] = "bone marrow"
            elif "PP005" in self.sample_fn or "PP006" in self.sample_fn:
                adata.obs["donor"] = "Donor1"
                adata.obs["organ"] = "lymph Node"
            elif "PP009" in self.sample_fn or "PP010" in self.sample_fn:
                adata.obs["donor"] = "Donor2"
                adata.obs["organ"] = "lung"
            elif "PP011" in self.sample_fn or "PP012" in self.sample_fn:
                adata.obs["donor"] = "Donor2"
                adata.obs["organ"] = "bone marrow"
            elif "PP013" in self.sample_fn or "PP014" in self.sample_fn:
                adata.obs["donor"] = "Donor2"
                adata.obs["organ"] = "lymph Node"
            adata.obs.index = self.sample_fn.split("_")[1].split("s")[0] + "nskept." + adata.obs.index
        adata.obs["cell_ontology_class"] = "Unknown"
        df1 = pd.read_csv(fn[1], sep="\t", index_col=0, header=None)
        df2 = pd.read_csv(fn[2], sep="\t", index_col=0, header=None)
        for i in df1.index:
            adata.obs["cell_ontology_class"].loc[i] = df1.loc[i][1]
        for i in df2.index:
            adata.obs["cell_ontology_class"].loc[i] = df2.loc[i][1]
        adata.X = scipy.sparse.csc_matrix(adata.X)

        return adata
