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
    this class. This dataloader only provides the subset of the published sata which has been made available through the
    covid-19 Cell Atlas.

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
        self.id = "human_colon_2019_10x_james_001_10.1038/s41590-020-0602-z"
        self.download_website = "https://covid19.cog.sanger.ac.uk/james20.processed.h5ad"
        self.organ = "colon"
        self.sub_tissue = "colonic immune cells"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Activated CD4 T': 'Activated CD4 T',
                'B cell IgA Plasma': 'B cell IgA Plasma',
                'B cell IgG Plasma': 'B cell IgG Plasma',
                'B cell cycling': 'B cell cycling',
                'B cell memory': 'B cell memory',
                'CD8 T': 'CD8 T',
                'Follicular B cell': 'Follicular',
                'ILC': 'ILC',
                'LYVE1 Macrophage': 'LYVE1 Macrophage',
                'Lymphoid DC': 'Lymphoid DC',
                'Macrophage': 'Macrophage',
                'Mast': 'Mast cell',
                'Monocyte': 'Monocyte',
                'NK': 'NK',
                'Tcm': 'Tcm',
                'Tfh': 'Tfh',
                'Th1': 'Th1',
                'Th17': 'Th17',
                'Treg': 'Treg',
                'cDC1': 'DC1',
                'cDC2': 'DC2',
                'cycling DCs': 'cycling DCs',
                'cycling gd T': 'cycling gd T',
                'gd T': 'gd T',
                'pDC': 'pDC',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/colon/james20.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[ADATA_IDS.lab] = "Teichmann"
        self.adata.uns[ADATA_IDS.year] = 2020
        self.adata.uns[ADATA_IDS.doi] = "10.1038/s41590-020-0602-z"
        self.adata.uns[ADATA_IDS.protocol] = '10x'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs['cell_type']
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs[ADATA_IDS.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col='gene_ids', new_index=ADATA_IDS.gene_id_ensembl)
