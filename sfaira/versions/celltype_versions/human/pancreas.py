from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_PANCREAS_V0 = {
    "names": {
        'Endocrine cell': ['Alpha cell', 'Beta cell', 'Gamma cell', 'Delta cell', 'Epsilon cell', 'Unclassified endocrine cell'],
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanPancreas(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_PANCREAS_V0
        }
        super(CelltypeVersionsHumanPancreas, self).__init__(**kwargs)
