from .external import CelltypeVersionsBase

# Version 0
ONTOLOGIES_MOUSE_OVARY_V0 = {
    "names": {
        'luteal cell': ['small luteal cell', 'large luteal cell'],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseOvary(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_OVARY_V0
        }
        super(CelltypeVersionsMouseOvary, self).__init__(**kwargs)
