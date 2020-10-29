from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_EYE_V0 = {
    "names": {
        'BPs': ['Retinal bipolar neuron type A', 'Retinal bipolar neuron type B', 'Retinal bipolar neuron type C', 'Retinal bipolar neuron type D'],
        'Rods': ['Retinal rod cell type A', 'Retinal rod cell type B', 'Retinal rod cell type C',]
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanEye(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_EYE_V0
        }
        super(CelltypeVersionsHumanEye, self).__init__(**kwargs)