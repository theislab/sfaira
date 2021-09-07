import anndata
import tarfile
import gzip
import scipy.io
import os
import pandas as pd
from sfaira.data import DatasetBase

SAMPLE_FNS = [
    "GSM3308545_NOD_08w_A",
    "GSM3308547_NOD_08w_C",
    "GSM3308548_NOD_14w_A",
    "GSM3308549_NOD_14w_B",
    "GSM3308550_NOD_14w_C",
    "GSM3308551_NOD_16w_A",
    "GSM3308552_NOD_16w_B",
    "GSM3308553_NOD_16w_C"
]


class Dataset(DatasetBase):
    """
    TODO add disease
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117770/suppl/GSE117770_RAW.tar"
        self.download_url_meta = f"private,{self.sample_fn}_annotation.csv"

        self.author = "Thompson"
        self.doi_journal = "10.1016/j.cmet.2019.01.021"
        self.normalization = "raw"
        self.organ = "pancreas"
        self.organism = "mouse"
        self.assay_sc = "10x 3' v2"
        self.state_exact = "diabetic"
        self.year = 2019
        self.sample_source = "primary_tissue"

        self.gene_id_symbols_var_key = "names"
        self.gene_id_ensembl_var_key = "ensembl"
        self.cell_type_obs_key = "celltypes"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    with tarfile.open(os.path.join(data_dir, 'GSE117770_RAW.tar')) as tar:
        for member in tar.getmembers():
            if "_matrix.mtx.gz" in member.name and sample_fn in member.name:
                name = "_".join(member.name.split("_")[:-1])
                with gzip.open(tar.extractfile(member), "rb") as mm:
                    x = scipy.io.mmread(mm).T.tocsr()
                obs = pd.read_csv(tar.extractfile(name + "_barcodes.tsv.gz"), compression="gzip", header=None,
                                  sep="\t", index_col=0)
                obs.index.name = None
                var = pd.read_csv(tar.extractfile(name + "_genes.tsv.gz"), compression="gzip", header=None,
                                  sep="\t")
                var.columns = ["ensembl", "names"]
                var.index = var["ensembl"].values
                adata = anndata.AnnData(X=x, obs=obs, var=var)
    adata.var_names_make_unique()
    celltypes = pd.read_csv(os.path.join(data_dir, sample_fn + "_annotation.csv"), index_col=0)
    adata = adata[celltypes.index]
    adata.obs["celltypes"] = celltypes

    return adata
