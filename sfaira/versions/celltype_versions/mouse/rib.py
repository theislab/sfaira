from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_RIB_V0 = [
    ["B cell", "CL:0000236"],
    ["endothelial cell", "CL:0000115"],
    ["macrophage", "CL:0000235"],
    ["stromal cell", "CL:0000499"],
    ["proliferative cell", "nan"],
    ["cartilage cell", "nan"],
    ["erythroblast", "nan"],
    ["granulocyte", "nan"],
    ["muscle cell", "nan"],
    ["neuron", "nan"],
    ["neutrophil", "nan"],
    ["osteoblast", "nan"],
    ["osteoclast", "nan"],
    ["oligodendrocyte", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_RIB_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseRib(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_RIB_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_RIB_V0
        }
        super(CelltypeVersionsMouseRib, self).__init__(**kwargs)
