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
        self.id = "human_lung_2020_10x_lukassen_002_10.1101/2020.03.13.991455"
        self.download_website = "https://covid19.cog.sanger.ac.uk/lukassen20_airway_orig.processed.h5ad"
        self.organ = "lung"
        self.sub_tissue = "bronchial epithelial cells"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Secretory3': 'Secretory',
                'Ciliated1': 'Multiciliated lineage',
                'Goblet': 'Secretory',
                'Ciliated2': 'Multiciliated lineage',
                'Club': 'Secretory',
                'Secretory2': 'Secretory',
                'FOXN4': 'Rare',
                'Basal1': 'Basal',
                'Secretory1': 'Secretory',
                'Fibroblast': '2_Fibroblast lineage',
                'Ionocyte': 'Rare',
                'Basal3': 'Basal',
                'Basal_Mitotic': 'Basal',
                'Basal2': 'Basal',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/lung/lukassen20_airway_orig.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['nCount_RNA'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[ADATA_IDS.lab] = 'Eils'
        self.adata.uns[ADATA_IDS.year] = 2020
        self.adata.uns[ADATA_IDS.doi] = "10.1101/2020.03.13.991455"
        self.adata.uns[ADATA_IDS.protocol] = '10x'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs['CellType']
        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs['state_exact'] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS.gene_id_ensembl)
