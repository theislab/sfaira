from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_ADIPOSE_V0 = {
    "names": {
        "lymphocyte": [
            "B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", 
            "myeloid cell", "NK cell"
        ],
        "T cell": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseAdipose(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_ADIPOSE_V0
        }
        super(CelltypeVersionsMouseAdipose, self).__init__(**kwargs)
