from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_SPLEEN_V0 = {
    "names": {
        'B cell (Plasmocyte)': ['Plasma_IgG', 'Plasma_IgM', 'Plasmablast'],
        'B cell': ['B_Hypermutation', 'B_follicular', 'B_mantle', 'B_T_doublet'],
        'Dendritic cell': ['DC_1', 'DC_2', 'DC_activated', 'DC_plasmacytoid'],
        'T cell': ["T_CD4_conv", "T_CD4_fh", "T_CD4_naive", "T_CD4_reg", "T_CD8_CTL", "T_CD8_MAIT", "T_CD8_activated",
                   "T_CD8_gd", ]
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanSpleen(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_SPLEEN_V0
        }
        super(CelltypeVersionsHumanSpleen, self).__init__(**kwargs)
