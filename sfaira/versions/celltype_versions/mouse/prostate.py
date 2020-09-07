from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_PROSTATE_V0 = [
    ["epithelial cell", "nan"],
    ["glandular epithelial cell", "nan"],
    ["T cell", "CL:0000084"],
    ["glandular cell", "nan"],
    ["stromal cell", "nan"],
    ["dendritic cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_PROSTATE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseProstate(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_PROSTATE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_PROSTATE_V0
        }
        super(CelltypeVersionsMouseProstate, self).__init__(**kwargs)
