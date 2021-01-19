from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_BLADDER_V0 = {
    "names": {
        "bladder cell": ["basal epithelial cell", "epithelial cell", "mesenchymal stromal cell", "smooth muscle cell",
                         "stromal cell", "umbrella cell"],
        "leukocyte": ["dendritic cell", "macrophage", "NK cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseBladder(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        print(__file__)
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_BLADDER_V0
        }
        super(CelltypeVersionsMouseBladder, self).__init__(**kwargs)
