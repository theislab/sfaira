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
        self.id = "human_rectum_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'rectum'
        self.sub_tissue = 'AdultRectum'
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
                'B cell': 'B cell',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'Dendritic cell': 'Dendritic cell',
                'Endothelial cell (APC)': 'Endothelial cell (APC)',
                'Enterocyte': 'Enterocyte',
                'Enterocyte progenitor': 'Enterocyte progenitor',
                'Epithelial cell': 'Epithelial cell',
                'Erythroid cell': 'Erythroid cell',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Macrophage': 'Macrophage',
                'Mast cell': 'Mast cell',
                'Monocyte': 'Monocyte',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Stromal cell': 'Stromal cell',
                'T cell': 'T cell',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "rectum", "hcl_AdultRectum_1.h5ad")
            self.adata = anndata.read(fn)

        self._convert_and_set_var_names(symbol_col="names", ensembl_col="ensembl")