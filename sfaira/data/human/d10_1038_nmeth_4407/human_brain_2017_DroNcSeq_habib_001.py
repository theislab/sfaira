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
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_brain_2017_DroNcSeq_habib_001_10.1038/nmeth.4407"
        self.download = "https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad"
        self.download_meta = None
        self.organ = "brain"
        self.sub_tissue = "hippocampus, prefrontal cortex"
        self.author = 'Regev'
        self.year = 2017
        self.doi = "10.1038/nmeth.4407"
        self.protocol = 'DroNcSeq'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'CellType'

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
        if fn is None:
            fn = os.path.join(self.path, "human", "brain", "habib17.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                   .multiply(1/10000)
