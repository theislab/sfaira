from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_RIB_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanRib(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_RIB_V0
        }
        super(CelltypeVersionsHumanRib, self).__init__(**kwargs)
