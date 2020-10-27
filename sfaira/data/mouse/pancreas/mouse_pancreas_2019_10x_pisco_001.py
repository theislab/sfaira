import anndata
import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA


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
        self.id = "mouse_pancreas_2019_10x_pisco_001_10.1101/661728"
        self.source = source
        if self.source == "aws":
            self.download_website = "https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/"
        elif self.source == "figshare":
            self.download_website = "https://ndownloader.figshare.com/articles/8273102/versions/2"
        else:
            raise ValueError("source %s not recognized" % self.source)
        self.organ = "pancreas"
        self.sub_tissue = "pancreas"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                "pancreatic ductal cel": "pancreatic ductal cell"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse/pancreas/tabula-muris-senis-droplet-processed-official-annotations-Pancreas.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse/pancreas/Pancreas_droplet.h5ad")
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

        self.adata.uns[ADATA_IDS.author] = "Quake"
        self.adata.uns[ADATA_IDS.year] = "2019"
        self.adata.uns[ADATA_IDS.doi] = "10.1101/661728"
        self.adata.uns[ADATA_IDS.protocol] = "10x"
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "mouse"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'norm'
        # self.adata.obs[ADATA_IDS.cell_ontology_class] is already set
        self.adata.obs[ADATA_IDS.cell_types_original] = self.adata.obs[ADATA_IDS.cell_ontology_class].values.tolist()
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs[ADATA_IDS.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS.gene_id_ensembl)
