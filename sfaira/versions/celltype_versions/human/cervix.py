from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_CERVIX_V0 = [
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fibroblast', "nan"],
    ['Loop of Henle', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_CERVIX_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanCervix(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_CERVIX_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_CERVIX_V0
        }
        super(CelltypeVersionsHumanCervix, self).__init__(**kwargs)
