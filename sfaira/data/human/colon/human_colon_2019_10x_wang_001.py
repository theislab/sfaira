import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS
import anndata
import numpy as np
import scipy.sparse


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_colon_2019_10x_wang_001_10.1084/jem.20191130"
        self.download_website = "https://covid19.cog.sanger.ac.uk/wang20_colon.processed.h5ad"
        self.organ = "colon"
        self.sub_tissue = "colon"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Progenitor': 'Enterocyte Progenitors',
                'Enterocyte': 'Enterocytes',
                'Goblet': 'Goblet cells',
                'TA': 'TA',
                'Paneth-like': 'Paneth cells',
                'Stem Cell': 'Stem cells',
                'Enteriendocrine': 'Enteroendocrine cells',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/colon/wang20_colon.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[ADATA_IDS.lab] = "Chen"
        self.adata.uns[ADATA_IDS.year] = 2019
        self.adata.uns[ADATA_IDS.doi] = "10.1084/jem.20191130"
        self.adata.uns[ADATA_IDS.protocol] = '10x'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs['CellType']
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs[ADATA_IDS.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS.gene_id_ensembl)
