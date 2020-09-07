from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_THYROID_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fasciculata cell', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Loop of Henle', "nan"],
    ['M2 Macrophage', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Primordial germ cell', "nan"],
    ['Proliferating T cell', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['Thyroid follicular cell', "nan"]
]

ONTOLOGIES_HUMAN_THYROID_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanThyroid(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_THYROID_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_THYROID_V0
        }
        super(CelltypeVersionsHumanThyroid, self).__init__(**kwargs)
