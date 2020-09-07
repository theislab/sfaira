from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_CHORIONICVILLUS_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['CB CD34+', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epithelial cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Loop of Henle', "nan"],
    ['M2 Macrophage', "nan"],
    ['Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"]
]
ONTOLOGIES_HUMAN_CHORIONICVILLUS_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanChorionicvillus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_CHORIONICVILLUS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_CHORIONICVILLUS_V0
        }
        super(CelltypeVersionsHumanChorionicvillus, self).__init__(**kwargs)
