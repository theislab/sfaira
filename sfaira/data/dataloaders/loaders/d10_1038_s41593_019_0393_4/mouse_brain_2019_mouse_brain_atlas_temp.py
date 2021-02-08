import anndata
import numpy as np
import os
import pandas
import zipfile
import scipy.io
from typing import Union
from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "mouse_brain_2019_10x_hove_001_10.1038/s41593-019-0393-4"

        self.download_url_data = \
            "https://www.brainimmuneatlas.org/data_files/toDownload/filtered_gene_bc_matrices_mex_WT_fullAggr.zip"
        self.download_url_meta = \
            "https://www.brainimmuneatlas.org/data_files/toDownload/annot_fullAggr.csv"

        self.author = "Movahedi"
        self.doi = "10.1038/s41593-019-0393-4"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "brain"
        self.organism = "mouse"
        self.protocol = "microwell-seq"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_ensembl_col = "ensembl"
        self.var_symbol_col = "names"

        self.obs_key_cellontology_class = self._ADATA_IDS_SFAIRA.cell_ontology_class
        self.obs_key_cellontology_id = self._ADATA_IDS_SFAIRA.cell_ontology_id
        self.obs_key_cellontology_original = self._ADATA_IDS_SFAIRA.cell_ontology_class

        self.class_maps = {
            "0": {
                "Microglia": "microglial cell",
                "T/NKT cells": "CD8-positive, alpha-beta T cell",
                "Monocytes": "monocyte"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "mouse", "temp_mouse_brain_atlas", "filtered_gene_bc_matrices_mex_WT_fullAggr.zip")
        fn_meta = os.path.join(self.path, "mouse", "temp_mouse_brain_atlas", "annot_fullAggr.csv")

        archive = zipfile.ZipFile(fn)
        x = scipy.io.mmread(archive.open('filtered_gene_bc_matrices_mex/mm10/matrix.mtx')).T.tocsr()
        self.adata = anndata.AnnData(x)
        var = pandas.read_csv(archive.open('filtered_gene_bc_matrices_mex/mm10/genes.tsv'), sep="\t", header=None)
        var.columns = ["ensembl", "name"]
        obs_names = pandas.read_csv(archive.open('filtered_gene_bc_matrices_mex/mm10/barcodes.tsv'),
                                    sep="\t",
                                    header=None
                                    )[0].values
        assert len(obs_names) == self.adata.shape[0]
        assert var.shape[0] == self.adata.shape[1]
        obs = pandas.read_csv(self.path + fn_meta)

        # Match annotation to raw data.
        obs.index = obs["cell"].values
        obs = obs.loc[[i in obs_names for i in obs.index], :]
        idx_tokeep = np.where([i in obs.index for i in obs_names])[0]
        self.adata = self.adata[idx_tokeep, :]
        obs_names = obs_names[idx_tokeep]
        idx_map = np.array([obs.index.tolist().index(i) for i in obs_names])
        self.adata = self.adata[idx_map, :]
        obs_names = obs_names[idx_map]

        # Assign attributes
        self.adata.obs_names = obs_names
        self.adata.var = var
        self.adata.obs = obs
        assert np.all(self.adata.obs_names == self.adata.obs["cell"].values)
