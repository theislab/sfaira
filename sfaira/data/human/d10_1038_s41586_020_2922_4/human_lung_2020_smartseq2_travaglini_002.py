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
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_lung_2020_smartseq2_travaglini_002_10.1038/s41586-020-2922-4"
        self.download = "https://www.synapse.org/#!Synapse:syn21041850"
        self.download_meta = None
        self.organ = "lung"
        self.sub_tissue = "proximal, medial, distal, blood"
        self.author = 'Krasnow'
        self.year = 2020
        self.doi = "10.1038/s41586-020-2922-4"
        self.protocol = 'smartseq2'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'

        self.class_maps = {
            "0": {
                'Adventitial Fibroblast_P1': 'Fibroblasts',
                'Adventitial Fibroblast_P2': 'Fibroblasts',
                'Adventitial Fibroblast_P3': 'Fibroblasts',
                'Airway Smooth Muscle_P1': 'Airway smooth muscle',
                'Airway Smooth Muscle_P2': 'Airway smooth muscle',
                'Airway Smooth Muscle_P3': 'Airway smooth muscle',
                'Alveolar Epithelial Type 1_P1': 'AT1',
                'Alveolar Epithelial Type 1_P2': 'AT1',
                'Alveolar Epithelial Type 1_P3': 'AT1',
                'Alveolar Epithelial Type 2_P1': 'AT2',
                'Alveolar Epithelial Type 2_P2': 'AT2',
                'Alveolar Epithelial Type 2_P3': 'AT2',
                'Alveolar Fibroblast_P1': 'Fibroblasts',
                'Alveolar Fibroblast_P2': 'Fibroblasts',
                'Alveolar Fibroblast_P3': 'Fibroblasts',
                'Artery_P1': 'Arterial',
                'Artery_P2': 'Arterial',
                'Artery_P3': 'Arterial',
                'B_P1': 'B cell lineage',
                'B_P2': 'B cell lineage',
                'B_P3': 'B cell lineage',
                'Basal_P1': 'Basal',
                'Basal_P2': 'Basal',
                'Basal_P3': 'Basal',
                'Basophil/Mast 1_P1': 'Mast cells',
                'Basophil/Mast 1_P2': 'Mast cells',
                'Basophil/Mast 1_P3': 'Mast cells',
                'Bronchial Vessel 1_P1': 'Bronchial Vessel 1',
                'CD4+ Memory/Effector T_P1': 'T cell lineage',
                'CD4+ Naive T_P1': 'T cell lineage',
                'CD4+ Naive T_P2': 'T cell lineage',
                'CD8+ Memory/Effector T_P1': 'T cell lineage',
                'CD8+ Naive T_P1': 'T cell lineage',
                'CD8+ Naive T_P2': 'T cell lineage',
                'Capillary Aerocyte_P1': 'Capillary',
                'Capillary Aerocyte_P2': 'Capillary',
                'Capillary Aerocyte_P3': 'Capillary',
                'Capillary Intermediate 1_P2': 'Capillary Intermediate 1',
                'Capillary_P1': 'Capillary',
                'Capillary_P2': 'Capillary',
                'Capillary_P3': 'Capillary',
                'Ciliated_P1': 'Multiciliated lineage',
                'Ciliated_P2': 'Multiciliated lineage',
                'Ciliated_P3': 'Multiciliated lineage',
                'Classical Monocyte_P1': 'Monocytes',
                'Club_P1': 'Secretory',
                'Club_P2': 'Secretory',
                'Club_P3': 'Secretory',
                'Dendritic_P1': 'Dendritic cells',
                'Differentiating Basal_P3': 'Basal',
                'Fibromyocyte_P3': 'Fibromyocyte',
                'Goblet_P1': 'Secretory',
                'Goblet_P2': 'Secretory',
                'Goblet_P3': 'Secretory',
                'IGSF21+ Dendritic_P2': 'Macrophages',
                'IGSF21+ Dendritic_P3': 'Macrophages',
                'Intermediate Monocyte_P2': 'Monocytes',
                'Intermediate Monocyte_P3': 'Monocytes',
                'Ionocyte_P3': 'Rare',
                'Lipofibroblast_P1': 'Fibroblasts',
                'Lymphatic_P1': 'Lymphatic EC',
                'Lymphatic_P2': 'Lymphatic EC',
                'Lymphatic_P3': 'Lymphatic EC',
                'Macrophage_P2': 'Macrophages',
                'Macrophage_P3': 'Macrophages',
                'Myeloid Dendritic Type 2_P3': 'Dendritic cells',
                'Myofibroblast_P2': 'Myofibroblasts',
                'Myofibroblast_P3': 'Myofibroblasts',
                'Natural Killer T_P2': 'T cell lineage',
                'Natural Killer T_P3': 'T cell lineage',
                'Natural Killer_P1': 'Innate lymphoid cells',
                'Natural Killer_P2': 'Innate lymphoid cells',
                'Natural Killer_P3': 'Innate lymphoid cells',
                'Neuroendocrine_P1': 'Rare',
                'Neuroendocrine_P3': 'Rare',
                'Neutrophil_P1': 'Monocytes',
                'Neutrophil_P2': 'Monocytes',
                'Neutrophil_P3': 'Monocytes',
                'Nonclassical Monocyte_P1': 'Monocytes',
                'Nonclassical Monocyte_P2': 'Monocytes',
                'Pericyte_P1': 'Fibroblasts',
                'Pericyte_P2': 'Fibroblasts',
                'Pericyte_P3': 'Fibroblasts',
                'Plasma_P3': 'B cell lineage',
                'Plasmacytoid Dendritic_P1': 'Dendritic cells',
                'Plasmacytoid Dendritic_P2': 'Dendritic cells',
                'Plasmacytoid Dendritic_P3': 'Dendritic cells',
                'Proliferating NK/T_P2': 'Innate lymphoid cells',
                'Proliferating NK/T_P3': 'Innate lymphoid cells',
                'Signaling Alveolar Epithelial Type 2_P1': 'AT2',
                'Signaling Alveolar Epithelial Type 2_P3': 'AT2',
                'Vascular Smooth Muscle_P1': '2_Smooth Muscle',
                'Vascular Smooth Muscle_P2': '2_Smooth Muscle',
                'Vascular Smooth Muscle_P3': '2_Smooth Muscle',
                'Vein_P2': 'Venous',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "lung", "facs_normal_lung_blood_scanpy.20200205.RC4.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = scipy.sparse.csc_matrix(self.adata.X)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['nReads'].values[:, None])) \
            .multiply(1 / 1000000)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
