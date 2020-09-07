from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_PROSTATE_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['Basal cell', "nan"],
    ['Club', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Epithelial cell (intermediated)', "nan"],
    ['Fasciculata cell', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Goblet cell', "nan"],
    ['Hillock', "nan"],
    ['Leukocytes', "nan"],
    ['Luminal', "nan"],
    ['Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Primordial germ cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_PROSTATE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanProstate(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_PROSTATE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_PROSTATE_V0
        }
        super(CelltypeVersionsHumanProstate, self).__init__(**kwargs)
