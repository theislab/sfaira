from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_BLOOD_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal chondrocyte', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['M2 Macrophage', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Proliferating T cell', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_BLOOD_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanBlood(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_BLOOD_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_BLOOD_V0
        }
        super(CelltypeVersionsHumanBlood, self).__init__(**kwargs)
