from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_LIMBMUSCLE_V0 = {
    "names": {
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseLimbmuscle(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LIMBMUSCLE_V0
        }
        super(CelltypeVersionsMouseLimbmuscle, self).__init__(**kwargs)
