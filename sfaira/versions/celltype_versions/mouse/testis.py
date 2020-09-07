from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_TESTIS_V0 = [
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
ONTOLOGIES_MOUSE_TESTIS_V0 = {
    "names": {
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseTestis(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_TESTIS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_TESTIS_V0
        }
        super(CelltypeVersionsMouseTestis, self).__init__(**kwargs)
