import anndata
import os
from typing import Union
import numpy as np
import scipy.sparse

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_brain_2017_DroNcSeq_habib_001_10.1038/nmeth.4407"

        self.download_url_data = "https://covid19.cog.sanger.ac.uk/habib17.processed.h5ad"
        self.download_url_meta = None

        self.author = "Regev"
        self.doi = "10.1038/nmeth.4407"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "brain"
        self.organism = "human"
        self.protocol = "DroNc-seq"
        self.state_exact = "healthy"
        self.year = 2017

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "CellType"

        self.class_maps = {
            "0": {
                "exPFC1": "Glutamatergic neurons from the PFC 1",
                "exPFC2": "Glutamatergic neurons from the PFC 2",
                "exDG": "Granule neurons from the hip dentate gyrus region",
                "GABA1": "GABAergic interneurons 1",
                "GABA2": "GABAergic interneurons 2",
                "exCA3": "Pyramidal neurons from the hip CA region 1",
                "exCA1": "Pyramidal neurons from the hip CA region 2",
                "ODC1": "Oligodendrocytes",
                "ASC1": "Astrocytes 1",
                "OPC": "Oligodendrocyte precursors",
                "ASC2": "Astrocytes 2",
                "Unclassified": "Unknown",
                "MG": "Microglia",
                "NSC": "Neuronal stem cells",
                "END": "Endothelial cells",
            },
        }

    def _load(self):
        fn = os.path.join(self.data_dir, "habib17.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["n_counts"].values[:, None]))\
                                   .multiply(1 / 10000)
