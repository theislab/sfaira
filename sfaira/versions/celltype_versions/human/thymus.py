from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_THYMUS_V0 = {
    "names": {
        'B cell': ['B_memory', 'B_naive', 'B_pro/pre'],
        'Dendritic cell': ['DC1', 'DC2'],
        'T cell': ['alpha_beta_T(entry)', 'gamma_delta_T', 'Treg', 'CD4+T', 'CD4+Tmem', 'CD8+T', 'CD8+Tmem']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanThymus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_THYMUS_V0
        }
        super(CelltypeVersionsHumanThymus, self).__init__(**kwargs)
