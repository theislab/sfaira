from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_SPLEEN_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B_Hypermutation', "nan"],
    ['B_T_doublet', "nan"],
    ['B_follicular', "nan"],
    ['B_mantle', "nan"],
    ['CB CD34+', "nan"],
    ['CD34_progenitor', "nan"],
    ['DC_1', "nan"],
    ['DC_2', "nan"],
    ['DC_activated', "nan"],
    ['DC_plasmacytoid', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fibroblast', "nan"],
    ['ILC', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['NK_CD160pos', "nan"],
    ['NK_FCGR3Apos', "nan"],
    ['NK_dividing', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Plasma_IgG', "nan"],
    ['Plasma_IgM', "nan"],
    ['Plasmablast', "nan"],
    ['Platelet', "nan"],
    ['Proliferating T cell', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T_CD4_conv', "nan"],
    ['T_CD4_fh', "nan"],
    ['T_CD4_naive', "nan"],
    ['T_CD4_reg', "nan"],
    ['T_CD8_CTL', "nan"],
    ['T_CD8_MAIT', "nan"],
    ['T_CD8_activated', "nan"],
    ['T_CD8_gd', "nan"],
    ['unknown', "nan"]
]
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
            "0": CELLTYPES_HUMAN_SPLEEN_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_SPLEEN_V0
        }
        super(CelltypeVersionsHumanSpleen, self).__init__(**kwargs)