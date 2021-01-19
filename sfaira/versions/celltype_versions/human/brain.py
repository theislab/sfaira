from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_BRAIN_V0 = {
    "names": {
        'Astrocyte': ['Astrocytes 1', 'Astrocytes 2'],
        'Fetal Neuron': ['Glutamatergic neurons from the PFC 1', 'Glutamatergic neurons from the PFC 2',
                         'Granule neurons from the hip dentate gyrus region', 'GABAergic interneurons 1',
                         'GABAergic interneurons 2', 'Pyramidal neurons from the hip CA region 1', 'Pyramidal neurons from the hip CA region 2']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanBrain(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_BRAIN_V0
        }
        super(CelltypeVersionsHumanBrain, self).__init__(**kwargs)
