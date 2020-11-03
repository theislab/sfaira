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
        self.id = "human_lung_2019_10x_braga_002_10.1038/s41591-019-0468-5"
        self.download_website = "https://covid19.cog.sanger.ac.uk/" \
                                "vieira19_Bronchi_anonymised.processed.h5ad"
        self.organ = "lung"
        self.sub_tissue = "bronchi"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Ciliated 1': 'Multiciliated lineage',
                'Club': 'Secretory',
                'Ciliated 2': 'Multiciliated lineage',
                'Ionocytes': 'Rare',
                'Basal 2': 'Basal',
                'Goblet_1': 'Secretory',
                'Goblet 2': 'Secretory',
                'Basal 1': 'Basal',
                'Dendritic cells': 'Dendritic cells',
                'B cells': 'B cell lineage',
                'Luminal_Macrophages': 'Macrophages',
                'Neutrophils': 'Monocytes',
                'Endothelial': '1_Endothelial',
                'Smooth muscle': '2_Smooth Muscle',
                'T and NK': '2_Lymphoid',
                'Fibroblasts': '2_Fibroblast lineage',
                'Lymphatic': 'Lymphatic EC',
                'Mast cells': 'Mast cells',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/lung/vieira19_Bronchi_anonymised.processed.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = np.expm1(self.adata.X)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Teichmann'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1038/s41591-019-0468-5"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'norm'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['CellType']
        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs['state_exact'] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)
