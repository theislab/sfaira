import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetMca


class Dataset(DatasetMca):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.id = "mouse_brain_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "brain"
        self.sub_tissue = "brain"

        self.class_maps = {
            "0": {
                'Astroglial cell(Bergman glia)(Brain)': 'Bergmann glial cell',
                'Astrocyte_Atp1b2 high(Brain)': 'astrocyte',
                'Astrocyte_Mfe8 high(Brain)': 'astrocyte',
                'Astrocyte_Pla2g7 high(Brain)': 'astrocyte',
                'Granulocyte_Ngp high(Brain)': 'granulocyte',
                'Hypothalamic ependymal cell(Brain)': 'ependymal cell',
                'Macrophage_Klf2 high(Brain)': 'macrophage',
                'Macrophage_Lyz2 high(Brain)': 'macrophage',
                'Microglia(Brain)': 'microglial cell',
                'Myelinating oligodendrocyte(Brain)': 'oligodendrocyte',
                'Oligodendrocyte precursor cell(Brain)': 'oligodendrocyte precursor cell',
                'Neuron(Brain)': 'neuron',
                'Pan-GABAergic(Brain)': 'GABAergic cell',
                'Schwann cell(Brain)': 'schwann cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Brain2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        self._load_mca(fn=fn, fn_meta=fn_meta)
