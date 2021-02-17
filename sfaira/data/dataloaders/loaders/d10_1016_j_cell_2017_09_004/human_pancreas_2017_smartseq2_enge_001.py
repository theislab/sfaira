import os
from typing import Union
import tarfile
import gzip
from io import StringIO
import anndata as ad
import pandas as pd
import scipy.sparse

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
        self.id = "human_pancreas_2017_smartseq2_enge_001_10.1016/j.cell.2017.09.004"

        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/suppl/GSE81547_RAW.tar"
        self.download_url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/matrix/GSE81547_series_matrix.txt.gz"

        self.author = "Quake"
        self.doi = "10.1016/j.cell.2017.09.004"
        self.healthy = True
        self.normalization = "raw"
        self.protocol = "Smart-seq2"
        self.organ = "islet of Langerhans"
        self.organism = "human"
        self.state_exact = "healthy"
        self.year = 2017

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "celltype"

        self.class_maps = {
            "0": {
                "alpha": "Alpha cell",
                "acinar": "Acinar cell",
                "ductal": "Ductal cell",
                "beta": "Beta cell",
                "unsure": "Unknown",
                "delta": "Delta cell",
                "mesenchymal": "Mesenchymal Cell"
            },
        }

    def _load(self):
        fn = [
            os.path.join(self.data_dir, "GSE81547_RAW.tar"),
            os.path.join(self.data_dir, "GSE81547_series_matrix.txt.gz")
        ]
        dfs = []
        with tarfile.open(fn[0]) as tar:
            for member in tar.getmembers():
                d = pd.read_csv(tar.extractfile(member), compression="gzip", header=None, sep="\t", index_col=0,
                                names=[member.name.split("_")[0]])
                dfs.append(d)
        adata = ad.AnnData(pd.concat(dfs, axis=1).iloc[1:-6].T)
        adata.X = scipy.sparse.csc_matrix(adata.X)
        with gzip.open(fn[1]) as f:
            file_content = [i.decode("utf-8") for i in f.readlines()]
        inputstring = ""
        for line in file_content:
            if "ID_REF" in line:
                inputstring += line
            if "!Sample_title" in line:
                inputstring += line[1:]
            if "!Sample_characteristics_ch1\t\"inferred_cell_type: alpha" in line:
                inputstring += line[1:]
        data = StringIO(inputstring)
        d = pd.read_csv(data, sep="\t").T
        d.columns = d.iloc[0]
        d.drop("Sample_title", inplace=True)
        d = d.reset_index().set_index("ID_REF")
        d.columns.name = None
        d.index.name = None
        adata.obs["celltype"] = [d.loc[i]["Sample_characteristics_ch1"].split(": ")[1] for i in adata.obs.index]
        adata.obs["patient"] = ["_".join(d.loc[i]["index"].split("_")[:2]) for i in adata.obs.index]

        return adata
