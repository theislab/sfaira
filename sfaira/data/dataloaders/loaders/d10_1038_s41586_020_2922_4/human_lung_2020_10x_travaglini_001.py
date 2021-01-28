import anndata
import os
from typing import Union
import scipy.sparse
import numpy as np

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "droplet_normal_lung_blood_scanpy.20200205.RC4.h5ad",
    "facs_normal_lung_blood_scanpy.20200205.RC4.h5ad"
]


class Dataset(DatasetBaseGroupLoadingManyFiles):
    """
    This data loader directly processes the data file provided under the download link.
    To obtain the file, you need to create a free account at https://www.synapse.org.
    You can then use those login credentials to download the file with python using the synapse client,
    installable via `pip install synapseclient`:

    import synapseclient
    import shutil
    syn = synapseclient.Synapse()
    syn.login("synapse_username","password")
    syn21625095 = syn.get(entity="syn21625095")
    shutil.move(syn21625095.path, "droplet_normal_lung_blood_scanpy.20200205.RC4.h5ad")

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            sample_fn: str,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        protocol = "10x" if self.sample_fn.split("_")[0] == "droplet" else "smartseq2"
        self.id = f"human_lung_2020_{protocol}_travaglini_{str(SAMPLE_FNS.index(self.sample_fn)).zfill(3)}_" \
                  f"10.1038/s41586-020-2922-4"

        self.download_url_data = "synapse,droplet_normal_lung_blood_scanpy.20200205.RC4.h5ad"
        self.download_url_meta = None

        self.author = "Krasnow"
        self.doi = "10.1038/s41586-020-2922-4"
        self.healthy = True
        self.normalization = "raw"
        self.organ = "lung"
        self.organism = "human"
        self.protocol = protocol
        self.state_exact = "healthy"
        self.year = 2020

        self.var_symbol_col = "index"

        self.class_maps = {
            "0": {
                "Adventitial Fibroblast_P1": "Fibroblasts",
                "Adventitial Fibroblast_P2": "Fibroblasts",
                "Adventitial Fibroblast_P3": "Fibroblasts",
                "Airway Smooth Muscle_P1": "Airway smooth muscle",
                "Airway Smooth Muscle_P2": "Airway smooth muscle",
                "Airway Smooth Muscle_P3": "Airway smooth muscle",
                "Alveolar Epithelial Type 1_P1": "AT1",
                "Alveolar Epithelial Type 1_P2": "AT1",
                "Alveolar Epithelial Type 1_P3": "AT1",
                "Alveolar Epithelial Type 2_P1": "AT2",
                "Alveolar Epithelial Type 2_P2": "AT2",
                "Alveolar Epithelial Type 2_P3": "AT2",
                "Alveolar Fibroblast_P1": "Fibroblasts",
                "Alveolar Fibroblast_P2": "Fibroblasts",
                "Alveolar Fibroblast_P3": "Fibroblasts",
                "Artery_P1": "Arterial",
                "Artery_P2": "Arterial",
                "Artery_P3": "Arterial",
                "B_P1": "B cell lineage",
                "B_P2": "B cell lineage",
                "B_P3": "B cell lineage",
                "Basal_P1": "Basal",
                "Basal_P2": "Basal",
                "Basal_P3": "Basal",
                "Basophil/Mast 1_P1": "Mast cells",
                "Basophil/Mast 1_P2": "Mast cells",
                "Basophil/Mast 1_P3": "Mast cells",
                "Basophil/Mast 2_P3": "Mast cells",
                "Bronchial Vessel 1_P1": "Bronchial Vessel 1",
                "Bronchial Vessel 1_P3": "Bronchial Vessel 1",
                "Bronchial Vessel 2_P1": "Bronchial Vessel 2",
                "Bronchial Vessel 2_P3": "Bronchial Vessel 2",
                "CD4+ Memory/Effector T_P1": "T cell lineage",
                "CD4+ Memory/Effector T_P2": "T cell lineage",
                "CD4+ Memory/Effector T_P3": "T cell lineage",
                "CD4+ Naive T_P1": "T cell lineage",
                "CD4+ Naive T_P2": "T cell lineage",
                "CD4+ Naive T_P3": "T cell lineage",
                "CD8+ Memory/Effector T_P1": "T cell lineage",
                "CD8+ Memory/Effector T_P2": "T cell lineage",
                "CD8+ Memory/Effector T_P3": "T cell lineage",
                "CD8+ Naive T_P1": "T cell lineage",
                "CD8+ Naive T_P2": "T cell lineage",
                "CD8+ Naive T_P3": "T cell lineage",
                "Capillary Aerocyte_P1": "Capillary",
                "Capillary Aerocyte_P2": "Capillary",
                "Capillary Aerocyte_P3": "Capillary",
                "Capillary Intermediate 1_P2": "Capillary Intermediate 1",
                "Capillary Intermediate 2_P2": "Capillary Intermediate 2",
                "Capillary_P1": "Capillary",
                "Capillary_P2": "Capillary",
                "Capillary_P3": "Capillary",
                "Ciliated_P1": "Multiciliated lineage",
                "Ciliated_P2": "Multiciliated lineage",
                "Ciliated_P3": "Multiciliated lineage",
                "Classical Monocyte_P1": "Monocytes",
                "Classical Monocyte_P2": "Monocytes",
                "Classical Monocyte_P3": "Monocytes",
                "Club_P1": "Secretory",
                "Club_P2": "Secretory",
                "Club_P3": "Secretory",
                "Differentiating Basal_P1": "Basal",
                "Differentiating Basal_P3": "Basal",
                "EREG+ Dendritic_P1": "Macrophages",
                "EREG+ Dendritic_P2": "Macrophages",
                "Fibromyocyte_P3": "Fibromyocyte",
                "Goblet_P3": "Secretory",
                "IGSF21+ Dendritic_P1": "Macrophages",
                "IGSF21+ Dendritic_P2": "Macrophages",
                "IGSF21+ Dendritic_P3": "Macrophages",
                "Intermediate Monocyte_P2": "Monocytes",
                "Ionocyte_P3": "Rare",
                "Lipofibroblast_P1": "Fibroblasts",
                "Lymphatic_P1": "Lymphatic EC",
                "Lymphatic_P2": "Lymphatic EC",
                "Lymphatic_P3": "Lymphatic EC",
                "Macrophage_P1": "Macrophages",
                "Macrophage_P2": "Macrophages",
                "Macrophage_P3": "Macrophages",
                "Mesothelial_P1": "Mesothelium",
                "Mucous_P2": "Submucosal Secretory",
                "Mucous_P3": "Submucosal Secretory",
                "Myeloid Dendritic Type 1_P1": "Dendritic cells",
                "Myeloid Dendritic Type 1_P2": "Dendritic cells",
                "Myeloid Dendritic Type 1_P3": "Dendritic cells",
                "Myeloid Dendritic Type 2_P1": "Dendritic cells",
                "Myeloid Dendritic Type 2_P2": "Dendritic cells",
                "Myeloid Dendritic Type 2_P3": "Dendritic cells",
                "Myofibroblast_P1": "Myofibroblasts",
                "Myofibroblast_P2": "Myofibroblasts",
                "Myofibroblast_P3": "Myofibroblasts",
                "Natural Killer T_P2": "T cell lineage",
                "Natural Killer T_P3": "T cell lineage",
                "Natural Killer_P1": "Innate lymphoid cells",
                "Natural Killer_P2": "Innate lymphoid cells",
                "Natural Killer_P3": "Innate lymphoid cells",
                "Neuroendocrine_P3": "Rare",
                "Nonclassical Monocyte_P1": "Monocytes",
                "Nonclassical Monocyte_P2": "Monocytes",
                "Nonclassical Monocyte_P3": "Monocytes",
                "OLR1+ Classical Monocyte_P2": "Monocytes",
                "Pericyte_P1": "Fibroblasts",
                "Pericyte_P2": "Fibroblasts",
                "Pericyte_P3": "Fibroblasts",
                "Plasma_P1": "B cell lineage",
                "Plasma_P3": "B cell lineage",
                "Plasmacytoid Dendritic_P1": "Dendritic cells",
                "Plasmacytoid Dendritic_P2": "Dendritic cells",
                "Plasmacytoid Dendritic_P3": "Dendritic cells",
                "Platelet/Megakaryocyte_P1": "Megakaryocytes",
                "Platelet/Megakaryocyte_P3": "Megakaryocytes",
                "Proliferating Basal_P1": "Basal",
                "Proliferating Basal_P3": "Basal",
                "Proliferating Macrophage_P1": "Macrophages",
                "Proliferating Macrophage_P2": "Macrophages",
                "Proliferating Macrophage_P3": "Macrophages",
                "Proliferating NK/T_P2": "Innate lymphoid cells",
                "Proliferating NK/T_P3": "Innate lymphoid cells",
                "Proximal Basal_P3": "Basal",
                "Proximal Ciliated_P3": "Multiciliated lineage",
                "Serous_P3": "Submucosal Secretory",
                "Signaling Alveolar Epithelial Type 2_P3": "AT2",
                "TREM2+ Dendritic_P1": "Macrophages",
                "TREM2+ Dendritic_P3": "Macrophages",
                "Vascular Smooth Muscle_P2": "2_Smooth Muscle",
                "Vascular Smooth Muscle_P3": "2_Smooth Muscle",
                "Vein_P1": "Venous",
                "Vein_P2": "Venous",
                "Vein_P3": "Venous",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "lung", self.sample_fn)
        if self.sample_fn.split("_")[0] == "droplet":
            norm_const = 1000000
        else:
            norm_const = 10000
        self.adata = anndata.read(fn)
        self.adata.X = scipy.sparse.csc_matrix(self.adata.X)
        self.adata.X = np.expm1(self.adata.X)
        self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs["nUMI"].values[:, None])) \
            .multiply(1 / norm_const)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
