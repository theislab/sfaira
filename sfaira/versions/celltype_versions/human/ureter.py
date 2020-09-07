from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_URETER_V0 = [
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Epithelial cell (intermediated)', "nan"],
    ['Fibroblast', "nan"],
    ['Intermediated cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_URETER_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanUreter(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_URETER_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_URETER_V0
        }
        super(CelltypeVersionsHumanUreter, self).__init__(**kwargs)
