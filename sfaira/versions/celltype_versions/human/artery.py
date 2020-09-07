from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_ARTERY_V0 = [
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epithelial cell', "nan"],
    ['Fibroblast', "nan"],
    ['M2 Macrophage', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Mesothelial cell', "nan"],
    ['Monocyte', "nan"],
    ['Myeloid cell', "nan"],
    ['Neutrophil', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_ARTERY_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanArtery(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_ARTERY_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_ARTERY_V0
        }
        super(CelltypeVersionsHumanArtery, self).__init__(**kwargs)
