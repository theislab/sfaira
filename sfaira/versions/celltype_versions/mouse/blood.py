from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_BLOOD_V0 = [
    ["B cell", "CL:0000236"],
    ["macrophage", "CL:0000235"],
    ["T cell", "CL:0000084"],
    ["NK cell", "nan"],
    ["neutrophil", "nan"],
    ["monocyte", "nan"],
    ["erythroblast", "nan"],
    ["dendritic cell", "nan"],
    ["basophil", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_BLOOD_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseBlood(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_BLOOD_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_BLOOD_V0
        }
        super(CelltypeVersionsMouseBlood, self).__init__(**kwargs)
