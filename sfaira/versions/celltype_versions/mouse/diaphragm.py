from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_DIAPHRAGM_V0 = [
    ["B cell", "CL:0000236"],
    ["endothelial cell", "CL:0000115"],
    ["macrophage", "CL:0000235"],
    ["mesenchymal stem cell", "CL:0000134"],
    ["skeletal muscle satellite cell", "CL:0000594"],
    ["T cell", "CL:0000084"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_DIAPHRAGM_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseDiaphragm(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_DIAPHRAGM_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_DIAPHRAGM_V0
        }
        super(CelltypeVersionsMouseDiaphragm, self).__init__(**kwargs)
