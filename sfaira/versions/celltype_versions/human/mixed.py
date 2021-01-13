from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_MIXED_V0 = [
    ['1.CD4rest', "nan"],
    ['10.CD8EM/TRMact', "nan"],
    ['10.CD8TEMRAact', "nan"],
    ['11.CD8TEMRA', "nan"],
    ['2.CD4act1', "nan"],
    ['2.CD4rest2', "nan"],
    ['3.CD4act1', "nan"],
    ['3.CD4act2', "nan"],
    ['4.CD4act2', "nan"],
    ['4.CD4act3', "nan"],
    ['5.CD4TRMrest', "nan"],
    ['5.CD4act3', "nan"],
    ['6.CD4TRMact', "nan"],
    ['6.CD4Treg', "nan"],
    ['7.CD4Treg', "nan"],
    ['7.CD8EM/TRMrest', "nan"],
    ['8.CD8EM/TRMact', "nan"],
    ['8.CD8EM/TRMrest', "nan"],
    ['9.CD8TEMRArest', "nan"],
    ['9.CD8TRMrest', "nan"],
    ['Unknown', "nan"]
]

ONTOLOGIES_HUMAN_MIXED_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanMixed(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_MIXED_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_MIXED_V0
        }
        super(CelltypeVersionsHumanMixed, self).__init__(**kwargs)
