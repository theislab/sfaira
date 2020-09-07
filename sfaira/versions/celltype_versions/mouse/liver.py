from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_LIVER_V0 = [
    ["B cell", "CL:0000236"],
    ["dendritic cell", "nan"],
    ["CD4-positive, alpha-beta T cell", "nan"],
    ["CD8-positive, alpha-beta T cell", "nan"],
    ["duct epithelial cell", "nan"],
    ["erythroblast", "nan"],
    ["endothelial cell of hepatic sinusoid", "CL:1000398"],
    ["granulocyte", "nan"],
    ["hepatic stellate cell", "CL:0000632"],
    ["hepatocyte", "CL:0000182"],
    ["Kupffer cell", "CL:0000091"],
    ["macrophage", "nan"],
    ["neutrophil", "CL:0000775"],
    ["NK cell", "CL:0000623"],
    ["plasmacytoid dendritic cell", "CL:0000784"],
    ["stromal cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_LIVER_V0 = {
    "names": {
        "myeloid leukocyte": [
            "macrophage", "dendritic cell", "plasmacytoid dendritic cell", "erythroblast",
            "granulocyte", "neutrophil"
        ],
        "T cell": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"],
        "mature NK T cell": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "NK cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseLiver(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_LIVER_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LIVER_V0
        }
        super(CelltypeVersionsMouseLiver, self).__init__(**kwargs)
