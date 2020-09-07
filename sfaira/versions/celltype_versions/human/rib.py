from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_RIB_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal Neuron', "nan"],
    ['Fetal chondrocyte', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Kidney intercalated cell', "nan"],
    ['Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['T cell', "nan"],
    ['hESC', "nan"]
]
ONTOLOGIES_HUMAN_RIB_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanRib(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_RIB_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_RIB_V0
        }
        super(CelltypeVersionsHumanRib, self).__init__(**kwargs)
