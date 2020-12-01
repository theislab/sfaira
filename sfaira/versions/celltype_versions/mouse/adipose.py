from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_ADIPOSE_V0 = [
    ["B cell", "CL:0000236"],
    ["CD4-positive, alpha-beta T cell", "nan"],
    ["CD8-positive, alpha-beta T cell", "nan"],
    ["endothelial cell", "CL:0000115"],
    ["epithelial cell", "CL:0000066"],
    ["erythroblast", "nan"],
    ["macrophage", "nan"],
    ["mesenchymal stem cell of adipose", "CL:0002570"],
    ["myeloid cell", "CL:0000763"],
    ["NK cell", "CL:0000623"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_ADIPOSE_V0 = {
    "names": {
        "lymphocyte": [
            "B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", 
            "myeloid cell", "NK cell"
        ],
        "T cell": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseAdipose(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_ADIPOSE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_ADIPOSE_V0
        }
        super(CelltypeVersionsMouseAdipose, self).__init__(**kwargs)
