from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_ADIPOSE_V0 = [
    ['B cell (Plasmocyte)', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epithelial cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Fibroblast', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_ADIPOSE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanAdipose(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_ADIPOSE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_ADIPOSE_V0
        }
        super(CelltypeVersionsHumanAdipose, self).__init__(**kwargs)
