from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_MARROW_V0 = {
    "names": {
        "granulocyte": ["basophil", "neutrophil", "mast cell"],
        "mature alpha-beta T cell": ["CD4-positive, alpha-beta T cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseMarrow(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_MARROW_V0
        }
        super(CelltypeVersionsMouseMarrow, self).__init__(**kwargs)
