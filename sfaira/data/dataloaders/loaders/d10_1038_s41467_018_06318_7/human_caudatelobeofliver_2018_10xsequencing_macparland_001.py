import anndata
import os
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "private,GSE115469.csv.gz"
        self.download_url_meta = "private,GSE115469_labels.txt"

        self.assay_sc = "10x 3' v2"
        self.author = "MacParland"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41467-018-06318-7"
        self.normalization = "raw"
        self.organ = "caudate lobe of liver"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.year = 2018

        self.gene_id_symbols_var_key = "index"
        self.cell_type_obs_key = "celltype"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE115469.csv.gz"),
        os.path.join(data_dir, "GSE115469_labels.txt")
    ]
    adata = anndata.read_csv(fn[0]).T
    celltype_df = pd.read_csv(fn[1], sep="\t").set_index("CellName")
    adata.obs["celltype"] = [str(celltype_df.loc[i]["Cluster#"]) for i in adata.obs.index]

    return adata
