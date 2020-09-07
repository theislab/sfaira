from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_LIMBMUSCLE_V0 = [
    ["B cell", "CL:0000236"],
    ["dendritic cell", "nan"],
    ["endothelial cell", "CL:0000115"],
    ["erythroblast", "nan"],
    ["macrophage", "CL:0000235"],
    ["mesenchymal stem cell", "CL:0000134"],
    ["monocyte progenitor", "nan"],
    ["muscle cell", "nan"],
    ["neutrophil", "nan"],
    ["Schwann cell", "CL:0002573"],
    ["smooth muscle cell", "CL:0000192"],
    ["stromal cell", "nan"],
    ["skeletal muscle cell", "CL:0000192"],
    ["skeletal muscle satellite cell", "CL:0000594"],
    ["T cell", "CL:0000084"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_LIMBMUSCLE_V0 = {
    "names": {
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseLimbmuscle(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_LIMBMUSCLE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LIMBMUSCLE_V0
        }
        super(CelltypeVersionsMouseLimbmuscle, self).__init__(**kwargs)
