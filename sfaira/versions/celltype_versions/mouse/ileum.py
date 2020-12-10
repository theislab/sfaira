from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_ILEUM_V0 = [
    ["B cell", "CL:0000236"],
    ["macrophage", "CL:0000235"],
    ["T cell", "CL:0000084"],
    ["dendritic cell", "nan"],
    ["mast cell", "nan"],
    ["paneth cell", "nan"],
    ["stromal cell", "nan"],
    ["epithelial cell", "nan"],
    ["epithelial cell villi", "nan"],
    ["enteroendocrine cell", "nan"],
    ["erythroblast", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_ILEUM_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseIleum(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_ILEUM_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_ILEUM_V0
        }
        super(CelltypeVersionsMouseIleum, self).__init__(**kwargs)
