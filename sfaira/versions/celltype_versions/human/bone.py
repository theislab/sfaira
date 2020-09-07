from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_BONE_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_BONE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanBone(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_BONE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_BONE_V0
        }
        super(CelltypeVersionsHumanBone, self).__init__(**kwargs)
