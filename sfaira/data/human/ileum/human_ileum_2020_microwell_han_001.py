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
        self.id = "human_ileum_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'ileum'
        self.sub_tissue = 'AdultIleum'
        self.dev_stage = 'Adult'
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
                'B cell': 'B cells',
                'B cell (Plasmocyte)': 'Plasma Cells',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell': 'Endothelial cell',
                'Endothelial cell (APC)': 'Endothelial cell',
                'Endothelial cell (endothelial to mesenchymal transition)': 'Endothelial cell',
                'Enterocyte': 'Enterocytes',
                'Enterocyte progenitor': 'Enterocytes',
                'Epithelial cell': 'Epithelial cell',
                'Fetal Neuron': 'Fetal neuron',
                'Fetal enterocyte': 'Enterocytes',
                'Fetal epithelial progenitor': 'Progenitors',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Fetal neuron': 'Fetal neuron',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Fibroblast': 'Fibroblasts',
                'Hepatocyte/Endodermal cell': 'Hepatocyte/Endodermal cell',
                'M2 Macrophage': 'M2 Macrophage',
                'Macrophage': 'Macrophage',
                'Mast cell': 'Mast cells',
                'Monocyte': 'Monocyte',
                'Neutrophil (RPS high)': 'Neutrophil (RPS high)',
                'Proliferating T cell': 'T cells',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cells',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "ileum", "hcl_AdultIleum_2.h5ad")
            self.adata = anndata.read(fn)
