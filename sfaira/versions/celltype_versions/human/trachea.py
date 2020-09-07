from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_TRACHEA_V0 = [
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['Chondrocyte', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Fetal chondrocyte', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Goblet cell', "nan"],
    ['Loop of Henle', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['Thyroid follicular cell', "nan"]
]

ONTOLOGIES_HUMAN_TRACHEA_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanTrachea(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_TRACHEA_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_TRACHEA_V0
        }
        super(CelltypeVersionsHumanTrachea, self).__init__(**kwargs)
