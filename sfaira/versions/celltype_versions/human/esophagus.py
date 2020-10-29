from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_ESOPHAGUS_V0 = {
    "names": {
        "Mono_macro": ["Monocyte", "Macrophage"],
        "B cell": ['B_CD27neg', 'B_CD27pos'],
        "T cell": ["T_CD4", "T_CD8", "NK_T_CD8_Cytotoxic"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanEsophagus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_ESOPHAGUS_V0
        }
        super(CelltypeVersionsHumanEsophagus, self).__init__(**kwargs)
