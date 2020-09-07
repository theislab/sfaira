from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_UTERUS_V0 = [
    ['AT2 cell', "nan"],
    ['B cell', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Loop of Henle', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Myeloid cell', "nan"],
    ['Primordial germ cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_UTERUS_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanUterus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_UTERUS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_UTERUS_V0
        }
        super(CelltypeVersionsHumanUterus, self).__init__(**kwargs)
