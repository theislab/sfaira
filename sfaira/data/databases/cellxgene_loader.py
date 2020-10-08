import numpy as np
import os
from typing import Union
from .external import DatasetBase
import anndata

from .external import DatasetGroupBase


class Dataset(DatasetBase):
    """
    This is a dataloader for downloaded h5ad from cellxgene.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None],
            fn: str,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.fn = fn
        self.species = str(fn).split("_")[2]
        self.id = str(fn).split(".")[0]
        self.organ = str(fn).split("_")[3]
        self.sub_tissue = None
        self.download_website = None  # TODO
        self.has_celltypes = True

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        fn = os.path.join(self.path, self.fn)
        adata = anndata.read(fn)
        adata.X = adata.raw.X

        self.adata.uns["lab"] = adata.uns["contributors"]["name"]
        self.adata.uns["year"] = None
        self.adata.uns["doi"] = None  # TODO
        if len(np.unique(adata.obs["organism"].values)) > 1:
            raise Warning("found multiple assay in data set %s" % self.fn)
        self.adata.uns["protocol"] = adata.obs["assay"].values[0]
        # Select tissue: blood is handled as a separate tissue in .obs
        #if len(np.unique(adata.obs["tissue"].values)) > 1:
        #    raise Warning("found multiple tissue in data set %s" % self.fn)
        #self.adata.uns["organ"] = adata.obs["tissue"].values[0]
        self.adata.uns["organ"] = str(self.fn).split("_")[3]
        if len(np.unique(adata.obs["organism"].values)) > 1:
            raise Warning("found multiple organisms in data set %s" % self.fn)
        self.adata.uns["animal"] = adata.obs["organism"].values[0]
        self.adata.uns["id"] = self.id
        self.adata.uns["wget_download"] = self.download_website
        self.adata.uns["has_celltypes"] = self.has_celltypes
        self.adata.uns["counts"] = 'raw'

        self.adata.obs["subtissue"] = self.sub_tissue
        self.adata.obs["dev_stage"] = adata.obs["development_stage"].values
        self.adata.obs["sex"] = adata.obs["sex"].values
        self.adata.obs["ethnicity"] = adata.obs["ethnicity"].values
        self.adata.obs["healthy"] = adata.obs["disease"].values == "normal"
        self.adata.obs["state_exact"] = adata.obs["disease"].values

        self.adata.obs["cell_ontology_id"] = adata.obs["cell_type_ontology_term_id"].values.tolist()
        self.adata.obs["cell_ontology_class"] = adata.obs["cell_type"].values.tolist()
        self.adata.obs["cell_types_original"] = adata.obs["free_annotation"].values.tolist()

        self._convert_and_set_var_names(symbol_col='names', ensembl_col='ensembl', new_index='ensembl')

