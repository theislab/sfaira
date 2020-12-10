import anndata
import os
from typing import Union
from .external import DatasetBase
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
        self.id = "human_brain_2017_DroNcSeq_habib_001_10.1038/nmeth.4407"
        self.download_website = "https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad"
        self.download_website_meta = None
        self.organ = "brain"
        self.sub_tissue = "hippocampus, prefrontal cortex"
        self.annotated = True

        self.class_maps = {
            "0": {
                'exPFC1': 'Glutamatergic neurons from the PFC 1',
                'exPFC2': 'Glutamatergic neurons from the PFC 2',
                'exDG': 'Granule neurons from the hip dentate gyrus region',
                'GABA1': 'GABAergic interneurons 1',
                'GABA2': 'GABAergic interneurons 2',
                'exCA3': 'Pyramidal neurons from the hip CA region 1',
                'exCA1': 'Pyramidal neurons from the hip CA region 2',
                'ODC1': 'Oligodendrocytes',
                'ASC1': 'Astrocytes 1',
                'OPC': 'Oligodendrocyte precursors',
                'ASC2': 'Astrocytes 2',
                'Unclassified': 'Unknown',
                'MG': 'Microglia',
                'NSC': 'Neuronal stem cells',
                'END': 'Endothelial cells',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "brain", "habib17.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = "Regev"
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = 2017
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = "10.1038/nmeth.4407"
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = 'DroNcSeq'
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['CellType']
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
