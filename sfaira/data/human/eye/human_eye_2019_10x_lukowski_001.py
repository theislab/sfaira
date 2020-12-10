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
        self.id = "human_eye_2019_10x_lukowski_001_10.15252/embj.2018100811"
        self.download_website = "https://covid19.cog.sanger.ac.uk/lukowski19.processed.h5ad"
        self.download_website_meta = None
        self.organ = "eye"
        self.sub_tissue = "retina"
        self.annotated = True

        self.class_maps = {
            "0": {
                'Muller cell': 'Muller cell',
                'amacrine cell': 'Amacrine cell',
                'microglial cell': 'Microglia',
                'retinal bipolar neuron type A': 'Retinal bipolar neuron type A',
                'retinal bipolar neuron type B': 'Retinal bipolar neuron type B',
                'retinal bipolar neuron type C': 'Retinal bipolar neuron type C',
                'retinal bipolar neuron type D': 'Retinal bipolar neuron type D',
                'retinal cone cell': 'Retinal cone cell',
                'retinal ganglion cell': 'Retinal ganglion cell',
                'retinal rod cell type A': 'Retinal rod cell type A',
                'retinal rod cell type B': 'Retinal rod cell type B',
                'retinal rod cell type C': 'Retinal rod cell type C',
                'unannotated': 'Unknown',
                'unspecified': 'Unknown',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "eye", "lukowski19.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Wong'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = '10.15252/embj.2018100811'
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

        self._convert_and_set_var_names(symbol_col='index', ensembl_col='gene_ids')
