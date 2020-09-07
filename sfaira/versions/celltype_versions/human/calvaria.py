from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_CALVARIA_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal chondrocyte', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal neuron', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Kidney intercalated cell', "nan"],
    ['Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_CALVARIA_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanCalvaria(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_CALVARIA_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_CALVARIA_V0
        }
        super(CelltypeVersionsHumanCalvaria, self).__init__(**kwargs)
