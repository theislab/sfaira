from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_UTERUS_V0 = [
    ["B cell", "CL:0000236"],
    ["dendritic cell", "nan"],
    ["granulocyte", "nan"],
    ["macrophage", "CL:0000235"],
    ["NK cell", "CL:0000623"],
    ["stromal cell", "nan"],
    ["endothelial cell", "nan"],
    ["glandular epithelial cell", "nan"],
    ["keratinocyte", "nan"],
    ["monocyte", "nan"],
    ["muscle cell", "nan"],
    ["smooth muscle cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_UTERUS_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseUterus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_UTERUS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_UTERUS_V0
        }
        super(CelltypeVersionsMouseUterus, self).__init__(**kwargs)
