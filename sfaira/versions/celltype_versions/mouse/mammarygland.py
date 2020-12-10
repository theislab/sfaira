from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_MAMMARYGLAND_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseMammarygland(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_MAMMARYGLAND_V0
        }
        super(CelltypeVersionsMouseMammarygland, self).__init__(**kwargs)
