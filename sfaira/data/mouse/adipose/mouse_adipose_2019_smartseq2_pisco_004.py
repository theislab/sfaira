import anndata
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
        self.id = "mouse_adipose_2019_smartseq2_pisco_004_10.1101/661728"
        self.source = source
        if self.source == "aws":
            self.download_website = "https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/"
        elif self.source == "figshare":
            self.download_website = "https://ndownloader.figshare.com/articles/8273102/versions/2"
        else:
            raise ValueError("source %s not recognized" % self.source)
        self.organ = "adipose"
        self.sub_tissue = "adipose"
        self.has_celltypes = True

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "adipose", "tabula-muris-senis-facs-processed-official-annotations-SCAT.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "adipose", "scat_facs.h5ad")
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

        self.adata.uns["lab"] = "Quake"
        self.adata.uns["year"] = "2019"
        self.adata.uns["doi"] = "10.1101/661728"
        self.adata.uns["protocol"] = "smartseq2"
        self.adata.uns["organ"] = self.organ
        self.adata.uns["subtissue"] = self.sub_tissue
        self.adata.uns["animal"] = "mouse"
        self.adata.uns["id"] = self.id
        self.adata.uns["wget_download"] = self.download_website
        self.adata.uns["has_celltypes"] = self.has_celltypes
        self.adata.uns["counts"] = 'norm'
        # self.adata.obs["cell_ontology_class"] is already set
        self.adata.obs["cell_types_original"] = self.adata.obs["cell_ontology_class"].values.tolist()
        self.adata.obs["healthy"] = True
        self.adata.obs["state_exact"] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index='ensembl')
