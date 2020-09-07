from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_RECTUM_V0 = [
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Enterocyte', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Enteroendocrine', "nan"],
    ['Erythroid cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Goblet', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Paneth-like', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stem Cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['TA', "nan"]
]
ONTOLOGIES_HUMAN_RECTUM_V0 = {
    "names": {
        'Epithelial cell': ['Paneth-like', 'Enteroendocrine', 'Goblet']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanRectum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_RECTUM_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_RECTUM_V0
        }
        super(CelltypeVersionsHumanRectum, self).__init__(**kwargs)