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
        self.id = "human_malegonad_2018_10x_guo_001_10.1038/s41422-018-0099-2"
        self.download_website = "https://covid19.cog.sanger.ac.uk/guo18_donor.processed.h5ad"
        self.download_website_meta = None
        self.organ = "malegonad"
        self.sub_tissue = "testis"
        self.annotated = True

        self.class_maps = {
            "0": {
                'Elongated Spermatids': 'Elongated Spermatids',
                'Leydig cells': 'Leydig cells',
                'Early Primary Spermatocytes': 'Early Primary Spermatocytes',
                'Round Spermatids': 'Round Spermatids',
                'Endothelial cells': 'Endothelial cells',
                'Macrophages': 'Macrophages',
                'Myoid cells': 'Myoid cells',
                'Differentiating Spermatogonia': 'Differentiating Spermatogonia',
                'Late primary Spermatocytes': 'Late primary Spermatocytes',
                'Spermatogonial Stem cell': 'Spermatogonial Stem cell',
                'Sertoli cells': 'Sertoli cells',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "malegonad", "guo18_donor.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = "Cairns"
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2018
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1038/s41422-018-0099-2"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['CellType']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
