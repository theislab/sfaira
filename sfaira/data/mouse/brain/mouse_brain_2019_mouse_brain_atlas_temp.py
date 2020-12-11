import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetMca


class Dataset(DatasetMca):

    id: str

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetMca.__init__(self=self, path=path, meta_path=meta_path, **kwargs)

        self.author = "Movahedi"
        self.doi = "10.1038/s41593-019-0393-4"
        self.download = \
            "www.brainimmuneatlas.org/data_files/toDownload/filtered_gene_bc_matrices_mex_WT_fullAggr.zip"
        self.download_meta = \
            "www.brainimmuneatlas.org/data_files/toDownload/annot_fullAggr.csv"
        self.healthy = True
        self.id = "mouse_brain_2019_10x_hove_001_10.1038/s41593-019-0393-4"
        self.normalization = 'raw'
        self.organ = "brain"
        self.protocol = "microwell"
        self.state_exact = "healthy"
        self.sub_tissue = "brain"
        self.year = "2019"

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
            fn = os.path.join(self.path, "mouse", "temp_mouse_brain_atlas", "matrix.mtx")
        fn_barcodes = os.path.join(self.path, "mouse", "temp_mouse_brain_atlas", "barcodes.tsv")
        fn_var = os.path.join(self.path, "mouse", "temp_mouse_brain_atlas", "genes.tsv")
        fn_meta = os.path.join(self.path, "mouse", "temp_mouse_brain_atlas", "annot_fullAggr.csv")

        self.adata = anndata.read_mtx(fn)
        self.adata = anndata.AnnData(self.adata.X.T)
        var = pandas.read_csv(fn_var, sep="\t", header=None)
        var.columns = ["ensembl", "name"]
        obs_names = pandas.read_csv(fn_barcodes, sep="\t", header=None)[0].values
        assert len(obs_names) == self.adata.shape[0]
        assert var.shape[0] == self.adata.shape[1]
        obs = pandas.read_csv(self.path + fn_meta)

        # Match annotation to raw data.
        obs.index = obs["cell"].values
        obs = obs.loc[[x in obs_names for x in obs.index], :]
        idx_tokeep = np.where([x in obs.index for x in obs_names])[0]
        self.adata = self.adata[idx_tokeep, :]
        obs_names = obs_names[idx_tokeep]
        idx_map = np.array([obs.index.tolist().index(x) for x in obs_names])
        self.adata = self.adata[idx_map, :]
        obs_names = obs_names[idx_map]

        # Assign attributes
        self.adata.obs_names = obs_names
        self.adata.var = var
        self.adata.obs = obs
        assert np.all(self.adata.obs_names == self.adata.obs["cell"].values)

        self.set_unkown_class_id(ids=["nan"])
