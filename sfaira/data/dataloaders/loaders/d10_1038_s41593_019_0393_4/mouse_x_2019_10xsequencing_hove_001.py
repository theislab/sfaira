import anndata
import numpy as np
import os
import pandas
import zipfile
import scipy.io
from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = \
            "https://www.brainimmuneatlas.org/data_files/toDownload/filtered_gene_bc_matrices_mex_WT_fullAggr.zip"
        self.download_url_meta = \
            "https://www.brainimmuneatlas.org/data_files/toDownload/annot_fullAggr.csv"

        self.assay_sc = "10x 3' v2"
        self.author = "Hove"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41593-019-0393-4"
        self.layer_counts = "X"
        self.organism = "Mus musculus"
        self.primary_data = True
        self.sample_source = "primary_tissue"
        self.state_exact = "healthy"
        self.year = 2019

        self.bio_sample_obs_key = "sample"
        self.cell_type_obs_key = "cluster"
        self.organ_obs_key = "organ"

        self.feature_id_var_key = "ensembl"
        self.feature_symbol_var_key = "name"
        self.feature_type = "rna"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    sample_organ_dict = {
        "Choroid plexus": "choroid plexus",
        "Dura mater": "dura mater",
        "Enr. SDM": "brain meninx",
        "Whole brain": "brain",
    }
    fn = [
        os.path.join(data_dir, "filtered_gene_bc_matrices_mex_WT_fullAggr.zip"),
        os.path.join(data_dir, "annot_fullAggr.csv")
    ]

    with zipfile.ZipFile(fn[0]) as archive:
        x = scipy.io.mmread(archive.open('filtered_gene_bc_matrices_mex/mm10/matrix.mtx')).T.tocsr()
        adata = anndata.AnnData(x)
        var = pandas.read_csv(archive.open('filtered_gene_bc_matrices_mex/mm10/genes.tsv'), sep="\t", header=None)
        var.columns = ["ensembl", "name"]
        obs_names = pandas.read_csv(archive.open('filtered_gene_bc_matrices_mex/mm10/barcodes.tsv'),
                                    sep="\t",
                                    header=None
                                    )[0].values
    obs = pandas.read_csv(fn[1])
    obs.fillna("isnan", inplace=True)

    # Match annotation to raw data.
    obs.index = obs["cell"].values
    obs = obs.loc[[i in obs_names for i in obs.index], :]
    idx_tokeep = np.where([i in obs.index for i in obs_names])[0]
    adata = adata[idx_tokeep, :]
    obs_names = obs_names[idx_tokeep]
    idx_map = np.array([obs.index.tolist().index(i) for i in obs_names])
    adata = adata[idx_map, :]
    obs_names = obs_names[idx_map]
    obs["organ"] = [sample_organ_dict[x] for x in obs["sample"].values]

    # Assign attributes
    adata.obs_names = obs_names
    adata.var = var
    adata.obs = obs

    return adata
