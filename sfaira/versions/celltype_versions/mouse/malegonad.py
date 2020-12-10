from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_MALEGONAD_V0 = [
    ["macrophage", "CL:0000235"],
    ["leydig cell", "nan"],
    ["elongating spermatid", "nan"],
    ["erythroblast", "nan"],
    ["pre-sertoli cell", "nan"],
    ["sertoli cell", "nan"],
    ["preleptotene spermatogonia", "nan"],
    ["spermatogonia", "nan"],
    ["spermatocyte", "nan"],
    ["spermatid", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_MALEGONAD_V0 = {
    "names": {
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseMalegonad(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_MALEGONAD_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_MALEGONAD_V0
        }
        super(CelltypeVersionsMouseMalegonad, self).__init__(**kwargs)
