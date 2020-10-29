from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_SPLEEN_V0 = {
    "names": {
        "T cell": [
            "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", 
            "immature NKT cell", "mature NK T cell"
        ]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseSpleen(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_SPLEEN_V0
        }
        super(CelltypeVersionsMouseSpleen, self).__init__(**kwargs)
