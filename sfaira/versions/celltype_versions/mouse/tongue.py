from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_TONGUE_V0 = [
    ["basal cell of epidermis", "CL:0002187"],
    ["keratinocyte", "CL:0000312"],
    ["Langerhans cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_TONGUE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseTongue(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_TONGUE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_TONGUE_V0
        }
        super(CelltypeVersionsMouseTongue, self).__init__(**kwargs)
