import anndata
import numpy as np
import os
from typing import Union
from .external import DatasetBase


class Dataset(DatasetBase):

    id: str

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            source: str = "aws",
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "mouse"
        self.id = "mouse_brain_2019_smartseq2_pisco_002_10.1101/661728"
        self.source = source
        if self.source == "aws":
            self.download = "https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/"
        elif self.source == "figshare":
            self.download = "https://ndownloader.figshare.com/articles/8273102/versions/2"
        else:
            raise ValueError("source %s not recognized" % self.source)
        self.organ = "brain"
        self.sub_tissue = "brain"
        self.annotated = True

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "brain", "tabula-muris-senis-facs-processed-official-annotations-Brain_Myeloid.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "brain", "Brain_Myeloid_facs.h5ad")
            else:
                raise ValueError("source %s not recognized" % self.source)
        self.adata = anndata.read_h5ad(fn)
        if self.source == "aws":
            self.adata.X = self.adata.raw.X
            self.adata.var = self.adata.raw.var
            del self.adata.raw
            self.adata.obsm = {}
            self.adata.varm = {}
            self.adata.uns = {}

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = "Quake"
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = "2019"
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = "10.1101/661728"
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = "smartseq2"
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = self.species
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.download_meta] = self.download_meta
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'norm'
        # self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] is already set
        self.set_unkown_class_id(ids=["nan", "Il6 expressing cells"])
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class].values.tolist()
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
