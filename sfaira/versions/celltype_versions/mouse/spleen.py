from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_SPLEEN_V0 = [
    ["B cell", "CL:0000236"],
    ["CD4-positive, alpha-beta T cell", "nan"],
    ["CD8-positive, alpha-beta T cell", "nan"],
    ["dendritic cell", "nan"],
    ["erythroblast", "CL:0000765"],
    ["granulocyte", "CL:0000094"],
    ["immature NKT cell", "nan"],
    ["macrophage", "nan"],
    ["macrophage dendritic cell progenitor", "CL:0002009"],
    ["mature NK T cell", "nan"],
    ["megakaryocyte-erythroid progenitor cell", "nan"],
    ["monocyte", "nan"],
    ["neutrophil", "nan"],
    ["NK cell", "CL:0000623"],
    ["plasma cell", "nan"],
    ["proerythroblast", "CL:0000547"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_SPLEEN_V0 = {
    "names": {
        "T cell": [
            "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", 
            "immature NKT cell", "mature NK T cell"
        ]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseSpleen(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_SPLEEN_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_SPLEEN_V0
        }
        super(CelltypeVersionsMouseSpleen, self).__init__(**kwargs)
