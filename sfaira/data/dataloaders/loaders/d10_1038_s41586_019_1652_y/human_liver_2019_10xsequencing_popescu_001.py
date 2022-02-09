import anndata
import os

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "private,fetal_liver_alladata_.h5ad"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Popescu"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41586-019-1652-y"
        self.layer_counts = "X"
        self.organ = "liver"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2019

        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"
        self.cell_type_obs_key = "cell.labels"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "fetal_liver_alladata_.h5ad")
    adata = anndata.read(fn)

    return adata
