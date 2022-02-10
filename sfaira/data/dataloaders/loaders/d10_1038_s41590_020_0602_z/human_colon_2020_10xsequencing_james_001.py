import anndata
import os
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://covid19.cog.sanger.ac.uk/james20.processed.h5ad"
        self.download_url_meta = None

        self.assay_sc_obs_key = "assay"
        self.author = "James"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41590-020-0602-z"
        self.individual_obs_key = "donor"
        self.layer_counts = "X"
        self.organ = "colon"
        self.organism = "Homo sapiens"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2020

        self.feature_id_var_key = "gene_ids"
        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"
        self.cell_type_obs_key = "cell_type"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "james20.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)
    # Assay maps are described here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7212050/
    adata.obs["assay"] = [{
        '290b': "10x 3' transcription profiling",
        '298c': "10x 3' transcription profiling",
        '302c': "10x 3' transcription profiling",
        '390c': "10x 5' transcription profiling",
        '417c': "10x 5' transcription profiling",
    }[x] for x in adata.obs["donor"].values]
    return adata
