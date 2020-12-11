import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase


class Dataset(DatasetBase):

    id: str

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "mouse"
        self.id = "mouse_brain_2019_10x_hove_001_10.1038/s41593-019-0393-4"
        self.download = \
            "www.brainimmuneatlas.org/data_files/toDownload/filtered_gene_bc_matrices_mex_WT_fullAggr.zip"
        self.download_meta = \
            "www.brainimmuneatlas.org/data_files/toDownload/annot_fullAggr.csv"
        self.organ = "brain"
        self.sub_tissue = "brain"
        self.annotated = True

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
        self._convert_and_set_var_names(symbol_col="names", ensembl_col="ensembl")
        self.adata.obs = obs
        assert np.all(self.adata.obs_names == self.adata.obs["cell"].values)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = "Movahedi"
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = "2019"
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = "10.1038/s41593-019-0393-4"
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = "microwell"
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = self.species
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.download_meta] = self.download_meta
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'
        # self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] is already set
        self.set_unkown_class_id(ids=["nan"])
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class].values.tolist()
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = "healthy"
