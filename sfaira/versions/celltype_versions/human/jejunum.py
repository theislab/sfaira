from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_JEJUNUM_V0 = [
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Enterocyte', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Fetal endocrine cell', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fibroblast', "nan"],
    ['Hepatocyte/Endodermal cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_JEJUNUM_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanJejunum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_JEJUNUM_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_JEJUNUM_V0
        }
        super(CelltypeVersionsHumanJejunum, self).__init__(**kwargs)
