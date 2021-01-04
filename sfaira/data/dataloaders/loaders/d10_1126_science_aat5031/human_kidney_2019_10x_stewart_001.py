import anndata
import os
from typing import Union
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the two raw data files which can be obtained from the `download_website`
    attribute of this class.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.organism = "loaders"
        self.id = "human_kidney_2019_10x_stewart_001_10.1126/science.aat5031"
        self.download = [
            'https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Mature_Full_v2.1.h5ad',
            'https://cellgeni.cog.sanger.ac.uk/BenKidney_v2.1/Fetal_full.h5ad'
        ]
        self.download_meta = None
        self.organ = "kidney"
        self.sub_tissue = "renal medulla, renal pelvis, ureter, cortex of kidney"
        self.author = 'Clatworthy'
        self.year = 2019
        self.doi = '10.1126/science.aat5031'
        self.protocol = '10x'
        self.normalization = 'norm'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.var_ensembl_col = 'ID'
        self.obs_key_cellontology_original = 'celltype'

        self.class_maps = {
            "0": {
                'Ascending vasa recta endothelium': 'Endothelial Cells - AVR',
                'B cell': 'B cell',
                'CD4 T cell': 'CD4 T cell',
                'CD8 T cell': 'CD8 T cell',
                'CNT/PC - proximal UB': 'CNT/PC - proximal UB',
                'Cap mesenchyme': 'Cap mesenchyme',
                'Connecting tubule': 'Connecting tubule',
                'Descending vasa recta endothelium': 'Endothelial Cells - AEA & DVR',
                'Distal S shaped body': 'Distal S shaped body',
                'Distal renal vesicle': 'Distal renal vesicle',
                'Distinct proximal tubule 1': 'Distinct proximal tubule 1',
                'Distinct proximal tubule 2': 'Distinct proximal tubule 2',
                'Endothelium': 'Endothelial Cells (unassigned)',
                'Epithelial progenitor cell': 'Epithelial progenitor',
                'Erythroid': 'Erythroid',
                'Fibroblast': 'Fibroblast',
                'Fibroblast 1': 'Fibroblast',
                'Fibroblast 2': 'Fibroblast',
                'Glomerular endothelium': 'Endothelial Cells - glomerular capillaries',
                'Indistinct intercalated cell': 'Indistinct intercalated cell',
                'Innate like lymphocyte': 'Innate like lymphocyte',
                'Loop of Henle': 'Loop of Henle',
                'MNP-a/classical monocyte derived': 'MNP-a/classical monocyte derived',
                'MNP-b/non-classical monocyte derived': 'MNP-b/non-classical monocyte derived',
                'MNP-c/dendritic cell': 'MNP-c/dendritic cell',
                'MNP-d/Tissue macrophage': 'MNP-d/Tissue macrophage',
                'Macrophage 1': 'Macrophage',
                'Macrophage 2': 'Macrophage',
                'Mast cell': 'Mast cell',
                'Mast cells': 'Mast cell',
                'Medial S shaped body': 'Medial S shaped body',
                'Megakaryocyte': 'Megakaryocyte',
                'Monocyte': 'Monocyte',
                'Myofibroblast': 'Myofibroblast',
                'Myofibroblast 1': 'Myofibroblast',
                'Myofibroblast 2': 'Myofibroblast',
                'NK cell': 'NK cell',
                'NKT cell': 'NKT cell',
                'Neuron': 'Neuron',
                'Neutrophil': 'Neutrophil',
                'Pelvic epithelium': 'Pelvic epithelium',
                'Pelvic epithelium - distal UB': 'Pelvic epithelium - distal UB',
                'Peritubular capillary endothelium 1': 'Peritubular capillary endothelium 1',
                'Peritubular capillary endothelium 2': 'Peritubular capillary endothelium 2',
                'Plasmacytoid dendritic cell': 'Plasmacytoid dendritic cell',
                'Podocyte': 'Podocyte',
                'Principal cell': 'Principal cell',
                'Proliferating B cell': 'Proliferating B cell',
                'Proliferating NK cell': 'Proliferating NK cell',
                'Proliferating Proximal Tubule': 'Proliferating Proximal Tubule',
                'Proliferating cDC2': 'Proliferating cDC2',
                'Proliferating cap mesenchyme': 'Proliferating cap mesenchyme',
                'Proliferating distal renal vesicle': 'Proliferating distal renal vesicle',
                'Proliferating fibroblast': 'Proliferating fibroblast',
                'Proliferating macrophage': 'Proliferating macrophage',
                'Proliferating monocyte': 'Proliferating monocyte',
                'Proliferating myofibroblast': 'Proliferating myofibroblast',
                'Proliferating stroma progenitor': 'Proliferating stroma progenitor',
                'Proximal S shaped body': 'Proximal S shaped body',
                'Proximal UB': 'Proximal UB',
                'Proximal renal vesicle': 'Proximal renal vesicle',
                'Proximal tubule': 'Proximal tubule',
                'Stroma progenitor': 'Stroma progenitor',
                'Thick ascending limb of Loop of Henle': 'Thick ascending limb of Loop of Henle',
                'Transitional urothelium': 'Transitional urothelium',
                'Type A intercalated cell': 'Type A intercalated cell',
                'Type B intercalated cell': 'Collecting Duct - Intercalated Cells Type B',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'pDC': 'pDC',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = [
                os.path.join(self.path, "human", "kidney", "Mature_Full_v2.1.h5ad"),
                os.path.join(self.path, "human", "kidney", "Fetal_full.h5ad")
            ]
        adult = anndata.read(fn[0])
        fetal = anndata.read(fn[1])
        adult.obs['development'] = 'adult'
        fetal.obs['development'] = 'fetal'
        self.adata = adult.concatenate(fetal)
        self.adata.X = np.expm1(self.adata.X)
