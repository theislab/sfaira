from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_DUODENUM_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanDuodenum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_DUODENUM_V0
        }
        super(CelltypeVersionsHumanDuodenum, self).__init__(**kwargs)
