import os
from typing import Union
from .external import DatasetHcl


class Dataset(DatasetHcl):
    """
    This is a dataloader for a the Human Cell Landscape dataset (Han et al. 2020. doi: 10.1038/s41586-020-2157-4).

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
        self.id = "human_placenta_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = 'Placenta'
        self.sub_tissue = 'Placenta'
        self.dev_stage = 'Fetus'
        self.class_maps = {
            "0": {
                'Fibroblast': 'Fibroblast',
                'Macrophage': 'Macrophage',
                'Epithelial cell': 'Epithelial cell',
                'Erythroid cell': 'Erythroid cell',
                'Fetal stromal cell': 'Fetal stromal cell',
                'Stromal cell': 'Stromal cell',
                'Smooth muscle cell': 'Smooth muscle cell',
                'Endothelial cell': 'Endothelial cell',
                'T cell': 'T cell',
                'Monocyte': 'Monocyte',
                'Neutrophil': 'Neutrophil',
                'Intermediated cell': 'Intermediated cell',
                'Dendritic cell': 'Dendritic cell',
                'CB CD34+': 'CB CD34+',
                'Stratified epithelial cell': 'Stratified epithelial cell',
                'Fetal neuron': 'Fetal neuron',
                'B cell (Plasmocyte)': 'B cell (Plasmocyte)',
                'Endothelial cell (APC)': 'Endothelial cell (APC)',
                'B cell': 'B cell',
                'Epithelial cell (intermediated)': 'Epithelial cell (intermediated)',
                'hESC': 'hESC',
                'Basal cell': 'Basal cell',
                'Fetal mesenchymal progenitor': 'Fetal mesenchymal progenitor',
                'Endothelial cell (endothelial to mesenchymal transition)': 'Endothelial cell (endothelial to mesenchymal transition)',
                'Erythroid progenitor cell (RP high)': 'Erythroid progenitor cell (RP high)',
                'Fetal epithelial progenitor': 'Fetal epithelial progenitor',
                'Fetal fibroblast': 'Fetal fibroblast',
                'Fetal skeletal muscle cell': 'Fetal skeletal muscle cell',
                'M2 Macrophage': 'M2 Macrophage',
                'Myeloid cell': 'Myeloid cell',
                'Antigen presenting cell (RPS high)': 'Antigen presenting cell (RPS high)',
            },
        }

    def _load(self, fn=None):
        self._load_hcl(fn=fn, sample_id="Placenta_1")
