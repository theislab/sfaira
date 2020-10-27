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
        self.id = "human_colon_2019_10x_smilie_001_10.1016/j.cell.2019.06.029"
        self.download_website = "https://covid19.cog.sanger.ac.uk/smillie19_epi.processed.h5ad"
        self.organ = "colon"
        self.sub_tissue = "colonic epithelium"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Cycling TA': 'Cycling TA',
                'TA 1': 'TA 1',
                'TA 2': 'TA 2',
                'Immature Enterocytes 2': 'Immature Enterocytes 2',
                'Immature Enterocytes 1': 'Immature Enterocytes 1',
                'Enterocyte Progenitors': 'Enterocyte Progenitors',
                'Immature Goblet': 'Immature Goblet',
                'Enterocytes': 'Enterocytes',
                'Secretory TA': 'Secretory TA',
                'Best4+ Enterocytes': 'Best4+ Enterocytes',
                'CD8+ IELs': 'CD8+ IELs',
                'Goblet': 'Goblet cells',
                'Stem': 'Stem cells',
                'Tuft': 'Tuft',
                'Follicular': 'Follicular',
                'Enteroendocrine': 'Enteroendocrine cells',
                'Plasma': 'Plasma Cells',
                'CD4+ Memory': 'CD4+ Memory',
                'CD8+ LP': 'CD8+ LP',
                'CD69- Mast': 'CD69- Mast',
                'Macrophages': 'Macrophage',
                'GC': 'Glial cells',
                'Cycling B': 'B cell cycling',
                'CD4+ Activated Fos-hi': 'CD4+ T Activated Fos-hi',
                'CD4+ Activated Fos-lo': 'CD4+ T Activated Fos-lo',
                'NKs': 'NK',
                'Cycling T': 'Cycling T',
                'M cells': 'M cells',
                'CD69+ Mast': 'CD69+ Mast',
                'MT-hi': 'MT-hi',
                'CD8+ IL17+': 'CD8+ IL17+',
                'CD4+ PD1+': 'CD4+ PD1+',
                'DC2': 'DC2',
                'Treg': 'Treg',
                'ILCs': 'ILC',
                'DC1': 'DC1',
                'WNT2B+ Fos-lo 1': 'WNT2B+ Fos-lo 1',
                'WNT5B+ 2': 'WNT5B+ 2',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/colon/smillie19_epi.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = "Regev"
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1016/j.cell.2019.06.029"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.animal] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['CellType']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)
