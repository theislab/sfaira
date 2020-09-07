from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_THYMUS_V0 = [
    ["abT cell", "nan"],
    ["B cell", "nan"],
    ["dendritic cell", "nan"],
    ["DN1 thymocyte", "nan"],
    ["DN2 thymocyte", "nan"],
    ["DN3 thymocyte", "nan"],
    ["DN4 thymocyte", "nan"],
    ["double positive T cell", "nan"],
    ["endothelial cell", "CL:0000115"],
    ["epithelial cell of thymus", "CL:0002293"],
    ["fibroblast", "nan"],
    ["gdT cell", "nan"],
    ["macrophage", "nan"],
    ["professional antigen presenting cell", "CL:0000145"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_THYMUS_V0 = {
    "names": {
        'double negative T cell': ["DN1 thymocyte", "DN2 thymocyte", "DN3 thymocyte", "DN4 thymocyte"],
        "immature T cell": [
            "DN1 thymocyte", "DN2 thymocyte", "DN3 thymocyte", "DN4 thymocyte", "double positive T cell"
        ],
        "mature T cell": ["abT cell", "gdT cell"],
        'thymocyte': [
            "DN1 thymocyte", "DN2 thymocyte", "DN3 thymocyte", "DN4 thymocyte", "double positive T cell",
            "gdT cell", "abT cell"
        ],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseThymus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_THYMUS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_THYMUS_V0
        }
        super(CelltypeVersionsMouseThymus, self).__init__(**kwargs)
