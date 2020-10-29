from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_THYMUS_V0 = {
    "names": {
        'double negative T cell': ["DN1 thymocyte", "DN2 thymocyte", "DN3 thymocyte", "DN4 thymocyte"],
        "immature T cell": [
            "DN1 thymocyte", "DN2 thymocyte", "DN3 thymocyte", "DN4 thymocyte", "double positive T cell"
        ],
        "mature T cell": ["abT cell", "gdT cell"],
        'thymocyte': [
            "DN1 thymocyte", "DN2 thymocyte", "DN3 thymocyte", "DN4 thymocyte", "double positive T cell",
            "gdT cell", "abT cell"
        ],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseThymus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_THYMUS_V0
        }
        super(CelltypeVersionsMouseThymus, self).__init__(**kwargs)
