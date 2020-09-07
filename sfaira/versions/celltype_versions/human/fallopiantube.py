from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_FALLOPIANTUBE_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fibroblast', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_FALLOPIANTUBE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanFallopiantube(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_FALLOPIANTUBE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_FALLOPIANTUBE_V0
        }
        super(CelltypeVersionsHumanFallopiantube, self).__init__(**kwargs)
