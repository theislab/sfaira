from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_SKIN_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['Basal cell 1', "nan"],
    ['Basal cell 2', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Epithelial cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal Neuron', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Kidney intercalated cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Proliferating T cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['WNT1', "nan"],
    ['channel', "nan"],
    ['folicular', "nan"],
    ['granular', "nan"],
    ['hESC', "nan"],
    ['melanocyte', "nan"],
    ['mitotic', "nan"],
    ['spinous', "nan"]
]
ONTOLOGIES_HUMAN_SKIN_V0 = {
    "names": {
        'immune': ['B cell', 'T cell', 'Dendritic cell', 'Erythroid cell', 'Erythroid progenitor cell (RP high)', 'Macrophage',
                   'Mast cell', 'Monocyte', 'Neutrophil', 'Neutrophil (RPS high)', 'Proliferating T cell'],
        'Basal cell': ['Basal cell 1', 'Basal cell 2']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanSkin(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_SKIN_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_SKIN_V0
        }
        super(CelltypeVersionsHumanSkin, self).__init__(**kwargs)
