from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_SKIN_V0 = [
    ["basal cell of epidermis", "CL:0002187"],
    ["bulge keratinocyte", "nan"],
    ["epidermal cell", "CL:0000362"],
    ["fibroblast", "nan"],
    ["keratinocyte stem cell", "CL:0002337"],
    ["macrophage", "nan"],
    ["stem cell of epidermis", "CL:1000428"],
    ["T cell", "CL:0000084"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_SKIN_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseSkin(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_SKIN_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_SKIN_V0
        }
        super(CelltypeVersionsMouseSkin, self).__init__(**kwargs)
