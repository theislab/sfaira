from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_PANCREAS_V0 = {
    "names": {
        "leukocyte": [
            "b cell",
            "dendritic cell",
            "granulocyte",
            "macrophage",
            "t cell"
        ],
        "endocrine cell": [
            "pancreatic A cell",
            "pancreatic D cell",
            "pancreatic PP cell"
        ],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMousePancreas(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_PANCREAS_V0
        }
        super(CelltypeVersionsMousePancreas, self).__init__(**kwargs)
