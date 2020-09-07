from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_STOMACH_V0 = [
    ["antral mucous cell", "nan"],
    ["dendritic cell", "nan"],
    ["G cell", "nan"],
    ["gastric mucosal cell", "nan"],
    ["epithelial cell", "nan"],
    ["muscle cell", "nan"],
    ["macrophage", "CL:0000235"],
    ["parietal cell", "nan"],
    ["pit cell", "nan"],
    ["proliferative cell", "nan"],
    ["stomach cell", "nan"],
    ["tuft cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_STOMACH_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseStomach(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_STOMACH_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_STOMACH_V0
        }
        super(CelltypeVersionsMouseStomach, self).__init__(**kwargs)
