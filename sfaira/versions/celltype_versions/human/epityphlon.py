from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_EPITYPHLON_V0 = [
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Enterocyte', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Epithelial cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_EPITYPHLON_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanEpityphlon(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_EPITYPHLON_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_EPITYPHLON_V0
        }
        super(CelltypeVersionsHumanEpityphlon, self).__init__(**kwargs)
