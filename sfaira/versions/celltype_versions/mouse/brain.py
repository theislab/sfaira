from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_BRAIN_V0 = {
    "names": {
        "T cell": ["CD8-positive, alpha-beta T cell", "yd T cells", "mature NK T cell"],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseBrain(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_BRAIN_V0
        }
        super(CelltypeVersionsMouseBrain, self).__init__(**kwargs)
