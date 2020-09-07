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
        self.id = "mouse_prostate_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "prostate"
        self.sub_tissue = "prostate"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Dendritic cell(Prostate)': 'dendritic cell',
                'Epithelial cell(Prostate)': 'epithelial cell',
                'Glandular epithelium(Prostate)': 'glandular epithelial cell',
                'Prostate gland cell(Prostate)': 'glandular cell',
                'Stromal cell(Prostate)': 'stromal cell',
                'T cell(Prostate)': 'T cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse/temp_mouse_atlas/500more_dge/Prostate2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse/temp_mouse_atlas/MCA_CellAssignments.csv")

        celltypes = pandas.read_csv(fn_meta, index_col=1)
        celltypes = celltypes.drop(['Unnamed: 0'], axis=1)

        data = pandas.read_csv(fn, sep=' ', header=0)
        self.adata = anndata.AnnData(data.T)
        self.adata = self.adata[np.array([x in celltypes.index for x in self.adata.obs_names])].copy()
        self.adata.obs = celltypes.loc[self.adata.obs_names, :]

        self.adata.uns["lab"] = "Guo"
        self.adata.uns["year"] = "2018"
        self.adata.uns["doi"] = "10.1016/j.cell.2018.02.001"
        self.adata.uns["protocol"] = "microwell-seq"
        self.adata.uns["organ"] = self.organ
        self.adata.uns["subtissue"] = self.sub_tissue  # TODO
        self.adata.uns["animal"] = "mouse"
        self.adata.uns["id"] = self.id
        self.adata.uns["wget_download"] = self.download_website
        self.adata.uns["has_celltypes"] = self.has_celltypes
        self.adata.uns["counts"] = 'raw'
        self.adata.obs["cell_ontology_class"] = self.adata.obs["Annotation"].values.tolist()
        self.set_unkown_class_id(ids=[np.nan, "nan"])
        self.adata.obs["cell_types_original"] = self.adata.obs["Annotation"].values.tolist()
        self.adata.obs["healthy"] = True
        self.adata.obs["state_exact"] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index='ensembl')
