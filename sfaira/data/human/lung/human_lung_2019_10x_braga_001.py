import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata
import numpy as np


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
        self.id = "human_lung_2019_10x_braga_001_10.1038/s41591-019-0468-5"
        self.download_website = "https://covid19.cog.sanger.ac.uk/" \
                                "vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad"
        self.organ = "lung"
        self.sub_tissue = "alveoli, parenchyma"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Ciliated 2': 'Multiciliated lineage',
                'Luminal_Macrophages': 'Macrophages',
                'Basal 1': 'Basal',
                'Dendritic cells': 'Dendritic cells',
                'Endothelial': '1_Endothelial',
                'Lymphatic': 'Lymphatic EC',
                'Ciliated 1': 'Multiciliated lineage',
                'Smooth muscle': '2_Smooth Muscle',
                'Type_1_alveolar': 'AT1',
                'Neutrophils': 'Monocytes',
                'Club': 'Secretory',
                'Basal 2': 'Basal',
                'B cells': 'B cell lineage',
                'T and NK': '2_Lymphoid',
                'Mesothelium': 'Mesothelium',
                'Mast cells': 'Mast cells',
                'Fibroblasts': '2_Fibroblast lineage',
                'Type 2 alveolar': 'AT2',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/lung/vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Teichmann'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1038/s41591-019-0468-5"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.animal] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'norm'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['CellType']
        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs['state_exact'] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)
