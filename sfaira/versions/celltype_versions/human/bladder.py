from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_BLADDER_V0 = [
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epithelial cell', "nan"],
    ['Epithelial cell (intermediated)', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Goblet cell', "nan"],
    ['Intermediated cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_BLADDER_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanBladder(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_BLADDER_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_BLADDER_V0
        }
        super(CelltypeVersionsHumanBladder, self).__init__(**kwargs)
