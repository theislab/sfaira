from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_BLADDER_V0 = [
    ["basal epithelial cell", "nan"],
    ["bladder urothelial cell", "CL:1001428"],
    ["dendritic cell", "nan"],
    ["endothelial cell", "CL:0000115"],
    ["epithelial cell", "nan"],
    ["macrophage", "nan"],
    ["mesenchymal stromal cell", "nan"],
    ["NK cell", "nan"],
    ["smooth muscle cell", "nan"],
    ["stromal cell", "nan"],
    ["umbrella cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_BLADDER_V0 = {
    "names": {
        "bladder cell": ["basal epithelial cell", "epithelial cell", "mesenchymal stromal cell", "smooth muscle cell",
                         "stromal cell", "umbrella cell"],
        "leukocyte": ["dendritic cell", "macrophage", "NK cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseBladder(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_BLADDER_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_BLADDER_V0
        }
        super(CelltypeVersionsMouseBladder, self).__init__(**kwargs)
