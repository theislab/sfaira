from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_THYMUS_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B_memory', "nan"],
    ['B_naive', "nan"],
    ['B_plasma', "nan"],
    ['B_pro/pre', "nan"],
    ['CB CD34+', "nan"],
    ['CD4+T', "nan"],
    ['CD4+Tmem', "nan"],
    ['CD8+T', "nan"],
    ['CD8+Tmem', "nan"],
    ['CD8αα', "nan"],
    ['DC1', "nan"],
    ['DC2', "nan"],
    ['DN', "nan"],
    ['DP', "nan"],
    ['ETP', "nan"],
    ['Endo', "nan"],
    ['Epi_GCM2', "nan"],
    ['Ery', "nan"],
    ['Fb_1', "nan"],
    ['Fb_2', "nan"],
    ['Fb_cycling', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['ILC3', "nan"],
    ['Lymph', "nan"],
    ['Mac', "nan"],
    ['Mast', "nan"],
    ['Mgk', "nan"],
    ['Mono', "nan"],
    ['NK', "nan"],
    ['NKT', "nan"],
    ['NMP', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Proliferating T cell', "nan"],
    ['T(agonist)', "nan"],
    ['TEC(myo)', "nan"],
    ['TEC(neuro)', "nan"],
    ['Treg', "nan"],
    ['VSMC', "nan"],
    ['aDC', "nan"],
    ['alpha_beta_T(entry)', "nan"],
    ['cTEC', "nan"],
    ['gamma_delta_T', "nan"],
    ['mTEC(I)', "nan"],
    ['mTEC(II)', "nan"],
    ['mTEC(III)', "nan"],
    ['mTEC(IV)', "nan"],
    ['mcTEC', "nan"],
    ['pDC', "nan"]
]
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
            "0": CELLTYPES_HUMAN_THYMUS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_THYMUS_V0
        }
        super(CelltypeVersionsHumanThymus, self).__init__(**kwargs)
