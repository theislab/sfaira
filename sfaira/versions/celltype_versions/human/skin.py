from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_SKIN_V0 = {
    "names": {
        'immune': ['B cell', 'T cell', 'Dendritic cell', 'Erythroid cell', 'Erythroid progenitor cell (RP high)', 'Macrophage',
                   'Mast cell', 'Monocyte', 'Neutrophil', 'Neutrophil (RPS high)', 'Proliferating T cell'],
        'Basal cell': ['Basal cell 1', 'Basal cell 2']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanSkin(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_SKIN_V0
        }
        super(CelltypeVersionsHumanSkin, self).__init__(**kwargs)
