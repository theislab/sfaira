import anndata
import os
from typing import Union
from .external import DatasetBase
import scipy.sparse
import numpy as np


class Dataset(DatasetBase):
    """
    This data loader directly processes the data file provided by the authors. To obtain the file, you need to create a
    free account at https://www.synapse.org. You can then use those login credentials to download the file with python
    using the synapse client, installable via `pip install synapseclient`:

    import synapseclient
    import shutil
    syn = synapseclient.Synapse()
    syn.login('synapse_username','password')
    syn21625142 = syn.get(entity='syn21625142')
    shutil.move(syn21625142.path, 'facs_normal_lung_blood_scanpy.20200205.RC4.h5ad')

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
        self.id = "human_lung_2020_smartseq2_travaglini_002_10.1038/s41586-020-2922-4"
        self.download_website = "https://www.synapse.org/#!Synapse:syn21041850"
        self.download_website_meta = None
        self.organ = "lung"
        self.sub_tissue = "proximal, medial, distal, blood"
        self.annotated = True

        self.class_maps = {
            "0": {
                'Intermediate Monocyte': 'Monocytes',
                'Adventitial Fibroblast': 'Fibroblasts',
                'Myofibroblast': 'Myofibroblasts',
                'Fibromyocyte': 'Fibromyocyte',
                'Basal': 'Basal',
                'IGSF21+ Dendritic': 'Macrophages',
                'CD8+ Memory/Effector T': 'T cell lineage',
                'CD4+ Naive T': 'T cell lineage',
                'Myeloid Dendritic Type 2': 'Dendritic cells',
                'Neuroendocrine': 'Rare',
                'Ciliated': 'Multiciliated lineage',
                'Nonclassical Monocyte': 'Monocytes',
                'Plasmacytoid Dendritic': 'Dendritic cells',
                'Vein': 'Venous',
                'Basophil/Mast 1': 'Mast cells',
                'Natural Killer T': 'T cell lineage',
                'Ionocyte': 'Rare',
                'Bronchial Vessel 1': 'Bronchial Vessel 1',
                'Natural Killer': 'Innate lymphoid cells',
                'Capillary Aerocyte': 'Capillary',
                'Vascular Smooth Muscle': '2_Smooth Muscle',
                'Macrophage': 'Macrophages',
                'Pericyte': 'Fibroblasts',
                'CD4+ Memory/Effector T': 'T cell lineage',
                'B': 'B cell lineage',
                'Lymphatic': 'Lymphatic EC',
                'Signaling Alveolar Epithelial Type 2': 'AT2',
                'Alveolar Epithelial Type 1': 'AT1',
                'Plasma': 'B cell lineage',
                'Lipofibroblast': 'Fibroblasts',
                'Capillary Intermediate 1': 'Capillary Intermediate 1',
                'Capillary': 'Capillary',
                'Alveolar Fibroblast': 'Fibroblasts',
                'Classical Monocyte': 'Monocytes',
                'Goblet': 'Secretory',
                'Airway Smooth Muscle': 'Airway smooth muscle',
                'Club': 'Secretory',
                'Proliferating NK/T': 'Innate lymphoid cells',
                'Alveolar Epithelial Type 2': 'AT2',
                'Differentiating Basal': 'Basal',
                'CD8+ Naive T': 'T cell lineage',
                'Artery': 'Arterial',
                'Neutrophil': 'Monocytes',
                'Dendritic': 'Dendritic cells',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "lung", "facs_normal_lung_blood_scanpy.20200205.RC4.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = scipy.sparse.csc_matrix(self.adata.X)
            self.adata.X = np.expm1(self.adata.X)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['nReads'].values[:, None])) \
                .multiply(1 / 1000000)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = 'Krasnow'
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = 2020
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = "10.1038/s41586-020-2922-4"
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = 'smartseq2'
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = ["_".join(i.split('_')[:-1]) for i in self.adata.obs['free_annotation']]
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class].astype('category')
        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.uns[self._ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
