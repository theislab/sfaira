from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_MUSCLE_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal Neuron', "nan"],
    ['Fetal chondrocyte', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['M2 Macrophage', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Primordial germ cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['Ventricle cardiomyocyte', "nan"],
    ['hESC', "nan"]
]
ONTOLOGIES_HUMAN_MUSCLE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanMuscle(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_MUSCLE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_MUSCLE_V0
        }
        super(CelltypeVersionsHumanMuscle, self).__init__(**kwargs)
