import anndata
import os
from typing import Union
from .external import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            source: str = "aws",
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "mouse"
        self.id = "mouse_lung_2019_10x_pisco_001_10.1101/661728"
        self.source = source
        if self.source == "aws":
            self.download_website = "https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/"
        elif self.source == "figshare":
            self.download_website = "https://ndownloader.figshare.com/articles/8273102/versions/2"
        else:
            raise ValueError("source %s not recognized" % self.source)
        self.organ = "lung"
        self.sub_tissue = "lung"
        self.annotated = True

        self.class_maps = {
            "0": {
                "ciliated columnar cell of tracheobronchial tree": "ciliated cell",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "lung", "tabula-muris-senis-droplet-processed-official-annotations-Lung.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "lung", "Lung_droplet.h5ad")
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
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = "10x"
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = "mouse"
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'norm'
        # self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] is already set
        self.set_unkown_class_id(ids=["nan"])
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class].values.tolist()
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
