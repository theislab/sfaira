from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_FEMALEGONAD_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Epithelial cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fasciculata cell', "nan"],
    ['Fetal Neuron', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal neuron', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Immature sertoli cell (Pre-Sertoli cell)', "nan"],
    ['Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['hESC', "nan"]
]
ONTOLOGIES_HUMAN_FEMALEGONAD_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanFemalegonad(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_FEMALEGONAD_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_FEMALEGONAD_V0
        }
        super(CelltypeVersionsHumanFemalegonad, self).__init__(**kwargs)
