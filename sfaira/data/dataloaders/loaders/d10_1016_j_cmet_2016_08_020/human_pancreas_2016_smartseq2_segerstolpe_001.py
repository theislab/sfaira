import anndata
import os
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    TODO: move state exact to diesase
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.processed.1.zip"
        self.download_url_meta = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt"

        self.author = "Segerstolpe"
        self.doi_journal = "10.1016/j.cmet.2016.08.020"
        self.normalization = "raw"
        self.organ = "pancreas"
        self.organism = "human"
        self.assay_sc = "Smart-seq2"
        self.year = 2016
        self.sample_source = "primary_tissue"

        self.gene_id_symbols_var_key = "index"

        self.cell_type_obs_key = "Characteristics[cell type]"
        self.state_exact_obs_key = "Characteristics[disease]"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "E-MTAB-5061.processed.1.zip"),
        os.path.join(data_dir, "E-MTAB-5061.sdrf.txt")
    ]
    df = pd.read_csv(fn[0], sep="\t")
    df.index = df.index.get_level_values(0)
    df = df.drop("#samples", axis=1)
    df = df.T.iloc[:, :26178]
    adata = anndata.AnnData(df)
    adata.obs = pd.read_csv(fn[1], sep="\t").set_index("Source Name").loc[adata.obs.index]
    # filter observations which are not cells (empty wells, low quality cells etc.)
    adata = adata[adata.obs["Characteristics[cell type]"] != "not applicable"].copy()

    return adata
