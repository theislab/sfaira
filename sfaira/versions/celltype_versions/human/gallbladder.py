from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_GALLBLADDER_V0 = [
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epithelial cell', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Goblet cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Myeloid cell', "nan"],
    ['Neutrophil', "nan"],
    ['Primordial germ cell', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['hESC', "nan"]
]
ONTOLOGIES_HUMAN_GALLBLADDER_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanGallbladder(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_GALLBLADDER_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_GALLBLADDER_V0
        }
        super(CelltypeVersionsHumanGallbladder, self).__init__(**kwargs)
