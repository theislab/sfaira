import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
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
        self.id = "human_ileum_2019_10x_martin_001_10.1016/j.cell.2019.08.008"
        self.download_website = "https://covid19.cog.sanger.ac.uk/martin19.processed.h5ad"
        self.organ = "ileum"
        self.sub_tissue = "ileum"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'T cells': 'T cells',
                'Plasma cells': 'Plasma Cells',
                'B cells': 'B cells',
                'MNP': 'MNP',
                'ILC': 'ILC',
                'Enterocytes': 'Enterocytes',
                'Fibs': 'Fibroblasts',
                'CD36+ endothelium': 'CD36+ endothelium',
                'Progenitors': 'Progenitors',
                'Goblets': 'Goblet cells',
                'Glial cells': 'Glial cells',
                'Cycling': 'Cycling',
                'ACKR1+ endothelium': 'ACKR1+ endothelium',
                'Pericytes': 'Pericytes',
                'Lymphatics': 'Lymphatics',
                'Mast cells': 'Mast cells',
                'SM': 'Smooth muscle cell',
                'TA': 'TA',
                'Paneth cells': 'Paneth cells',
                'Enteroendocrines': 'Enteroendocrine cells',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/ileum/martin19.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)
            self.adata = self.adata[self.adata.obs['CellType'] != 'Doublets'].copy()

        self.adata.uns[ADATA_IDS_SFAIRA.author] = "Kenigsberg"
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1016/j.cell.2019.08.008"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['CellType']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col='gene_ids', new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)
