from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_LIVER_V0 = {
    "names": {
        "myeloid leukocyte": [
            "macrophage", "dendritic cell", "plasmacytoid dendritic cell", "erythroblast",
            "granulocyte", "neutrophil"
        ],
        "T cell": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"],
        "mature NK T cell": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "NK cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseLiver(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LIVER_V0
        }
        super(CelltypeVersionsMouseLiver, self).__init__(**kwargs)
