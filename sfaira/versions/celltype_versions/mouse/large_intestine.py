from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_LARGEINTESTINE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseLargeintestine(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LARGEINTESTINE_V0
        }
        super(CelltypeVersionsMouseLargeintestine, self).__init__(**kwargs)
