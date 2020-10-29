from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_RECTUM_V0 = {
    "names": {
        'Epithelial cell': ['Paneth-like', 'Enteroendocrine', 'Goblet']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanRectum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_RECTUM_V0
        }
        super(CelltypeVersionsHumanRectum, self).__init__(**kwargs)