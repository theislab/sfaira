import anndata
import os
from typing import Union
from .external import DatasetBase


class Dataset(DatasetBase):
    """
    This is a dataloader for a the Human Cell Landscape dataset (Han et al. 2020. doi: 10.1038/s41586-020-2157-4).
    In order to obtain the required preprocessed datafiles, please use the notebook provided in this repository under:
    sfaira/data/download_scripts/get_and_preprocess_HumanCellLandscape.ipynb

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
        self.id = "human_brain_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'brain'
        self.sub_tissue = 'FetalBrain'
        self.dev_stage = 'Fetus'
        self.download = 'https://figshare.com/articles/HCL_DGE_Data/7235471'
        self.download_meta = None
        self.author = 'Guo'
        self.year = 2020
        self.doi = '10.1038/s41586-020-2157-4'
        self.protocol = 'microwell'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'names'
        self.var_ensembl_col = 'ensembl'
        self.obs_key_cellontology_original = 'cell_ontology_class'

        self.class_maps = {
            "0": {
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
                'Astrocyte': 'Astrocyte',
                'B cell': 'B cell',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'CB CD34+': 'CB CD34+',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cells',
                'Endothelial cell (APC)': 'Endothelial cells',
                'Erythroid cell': 'Erythroid cell',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Fetal Neuron': 'Fetal Neuron',
                'Fetal endocrine cell': 'Fetal endocrine cell',
                'Fetal enterocyte ': 'Fetal enterocyte ',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Fetal neuron': 'Fetal Neuron',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Fibroblast': 'Fibroblast',
                'Gastric endocrine cell': 'Gastric endocrine cell',
                'Goblet cell': 'Goblet cell',
                'Macrophage': 'Macrophage',
                'Monocyte': 'Monocyte',
                'Neutrophil': 'Neutrophil',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Oligodendrocyte': 'Oligodendrocytes',
                'Primordial germ cell': 'Primordial germ cell',
                'Sinusoidal endothelial cell': 'Sinusoidal endothelial cell',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
                'hESC': 'Neuronal stem cells'
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "brain", "hcl_FetalBrain_4.h5ad")
            self.adata = anndata.read(fn)
