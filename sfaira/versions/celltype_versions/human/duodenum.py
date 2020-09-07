from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_DUODENUM_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Enterocyte', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Fetal endocrine cell', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fibroblast', "nan"],
    ['Goblet cell', "nan"],
    ['Hepatocyte/Endodermal cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_DUODENUM_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanDuodenum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_DUODENUM_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_DUODENUM_V0
        }
        super(CelltypeVersionsHumanDuodenum, self).__init__(**kwargs)
