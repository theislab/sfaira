import anndata
import numpy as np
import os
import pandas
import zipfile
import scipy.io


def load(data_dir, **kwargs):
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
    obs["organ"] = obs["sample"].values

    # Assign attributes
    adata.obs_names = obs_names
    adata.var = var
    adata.obs = obs

    return adata
