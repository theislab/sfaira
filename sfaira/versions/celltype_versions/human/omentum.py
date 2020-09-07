from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_OMENTUM_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epithelial cell', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fibroblast', "nan"],
    ['M2 Macrophage', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Mesothelial cell', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_OMENTUM_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanOmentum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_OMENTUM_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_OMENTUM_V0
        }
        super(CelltypeVersionsHumanOmentum, self).__init__(**kwargs)
