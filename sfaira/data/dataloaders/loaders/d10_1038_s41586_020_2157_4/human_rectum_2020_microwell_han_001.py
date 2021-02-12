from typing import Union
from .base import Dataset_d10_1038_s41586_020_2157_4


class Dataset(Dataset_d10_1038_s41586_020_2157_4):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.id = "human_rectum_2020_microwell_han_001_10.1038/s41586-020-2157-4"
        self.organ = "rectum"
        self.class_maps = {
            "0": {
                "B cell": "B cell",
                "B cell (Plasmocyte)": "B cell (Plasmocyte)",
                "Dendritic cell": "Dendritic cell",
                "Endothelial cell (APC)": "Endothelial cell (APC)",
                "Enterocyte": "Enterocyte",
                "Enterocyte progenitor": "Enterocyte progenitor",
                "Epithelial cell": "Epithelial cell",
                "Erythroid cell": "Erythroid cell",
                "Fetal stromal cell": "Fetal stromal cell",
                "Macrophage": "Macrophage",
                "Mast cell": "Mast cell",
                "Monocyte": "Monocyte",
                "Smooth muscle cell": "Smooth muscle cell",
                "Stromal cell": "Stromal cell",
                "T cell": "T cell",
            },
        }

    def _load(self):
        self._load_generalized(sample_id="AdultRectum_1")
