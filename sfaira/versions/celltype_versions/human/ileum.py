from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_ILEUM_V0 = [
    ['ACKR1+ endothelium', "nan"],
    ['B cells', "nan"],
    ['CD36+ endothelium', "nan"],
    ['Cycling', "nan"],
    ['Dendritic cell', "nan"],
    ['Enterocytes', "nan"],
    ['Enteroendocrine cells', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal neuron', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblasts', "nan"],
    ['Glial cells', "nan"],
    ['Goblet cells', "nan"],
    ['Hepatocyte/Endodermal cell', "nan"],
    ['ILC', "nan"],
    ['Lymphatics', "nan"],
    ['M2 Macrophage', "nan"],
    ['MNP', "nan"],
    ['Macrophage', "nan"],
    ['Mast cells', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Paneth cells', "nan"],
    ['Pericytes', "nan"],
    ['Plasma Cells', "nan"],
    ['Progenitors', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stem Cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cells', "nan"],
    ['TA', "nan"]
]
ONTOLOGIES_HUMAN_ILEUM_V0 = {
    "names": {
        'Endothelial cell': ['ACKR1+ endothelium', 'CD36+ endothelium'],
        'Epithelial cell': ['Goblet cells', 'Enterocytes', 'Paneth cells', 'Enteroendocrine cells']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanIleum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_ILEUM_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_ILEUM_V0
        }
        super(CelltypeVersionsHumanIleum, self).__init__(**kwargs)
