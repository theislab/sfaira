from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_ILEUM_V0 = {
    "names": {
        'Endothelial cell': ['ACKR1+ endothelium', 'CD36+ endothelium'],
        'Epithelial cell': ['Goblet cells', 'Enterocytes', 'Paneth cells', 'Enteroendocrine cells']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanIleum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_ILEUM_V0
        }
        super(CelltypeVersionsHumanIleum, self).__init__(**kwargs)
