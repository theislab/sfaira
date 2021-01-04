import anndata
import numpy as np
import os
import pandas
from typing import Union
from .base import Dataset_d10_1016_j_cell_2018_02_001


class Dataset(Dataset_d10_1016_j_cell_2018_02_001):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)

        self.id = "mouse_kidney_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "kidney"
        self.sub_tissue = "kidney"

        self.class_maps = {
            "0": {
                'Adipocyte(Fetal_Kidney)': 'fetal adipocyte',
                'B cell(Kidney)': 'B cell',
                'Dendritic cell_Ccr7 high(Kidney)': 'dendritic cell',
                'Dendritic cell_Cst3 high(Kidney)': 'dendritic cell',
                'Distal collecting duct principal cell_Cldn4 high(Kidney)': 'kidney collecting duct principal cell',
                'Distal collecting duct principal cell_Hsd11b2 high(Kidney)': 'kidney collecting duct principal cell',
                'Distal convoluted tubule_Pvalb high(Kidney)': 'kidney distal convoluted tubule epithelial cell',
                'Distal convoluted tubule_S100g high(Kidney)': 'kidney distal convoluted tubule epithelial cell',
                'Endothelial cell(Kidney)': 'fenestrated cell',
                'Epithelial cell_Cryab high(Kidney)': "epithelial cell",
                'Fenestrated endothelial cell_Plvap high(Kidney)': 'fenestrated cell',
                'Fenestrated endothelial cell_Tm4sf1 high(Kidney)': 'fenestrated cell',
                'Glomerular epithelial cell_Aldh1a2 high(Fetal_Kidney)': 'glomerular epithelial cell',
                'Intercalated cells of collecting duct_Aqp6 high(Kidney)': 'kidney collecting duct epithelial cell',
                'Intercalated cells of collecting duct_Slc26a4 high(Kidney)': 'kidney collecting duct epithelial cell',
                'Macrophage_Ccl4 high (Kidney)': 'macrophage',
                'Macrophage_Lyz2 high(Kidney)': 'macrophage',
                'Metanephric mesenchyme(Fetal_Kidney)': 'fetal mesenchymal cell',
                'Neutrophil progenitor_S100a8 high(Kidney)': 'neutrophil progenitor',
                'Proximal tubule brush border cell(Kidney)': 'brush cell',
                'Proximal tubule cell_Cyp4a14 high(Kidney)': 'epithelial cell of proximal tubule',
                'Proximal tubule cell_Osgin1 high(Kidney)': 'epithelial cell of proximal tubule',
                'S1 proximal tubule cells(Kidney)': 'epithelial cell of proximal tubule',
                'S3 proximal tubule cells(Kidney)': 'epithelial cell of proximal tubule',
                'Stromal cell_Ankrd1 high(Kidney)': 'fibroblast',
                'Stromal cell_Cxcl10 high(Kidney)': 'fibroblast',
                'Stromal cell_Dcn high(Kidney)': 'fibroblast',
                'Stromal cell_Mgp high(Fetal_Kidney)': 'fibroblast',
                'Stromal cell_Mgp high(Kidney)': 'fibroblast',
                'Stromal cell_Ptgds high(Kidney)': 'fibroblast',
                'T cell(Kidney)': 'T cell',
                'Thick ascending limb of the loop of Henle(Kidney)': 'kidney loop of Henle ascending limb epithelial cell',
                'Ureteric epithelium(Kidney)': 'ureteric epithelial cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Kidney2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_generalized(fn=fn, fn_meta=fn_meta)
