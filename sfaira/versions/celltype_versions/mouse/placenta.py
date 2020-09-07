from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_PLACENTA_V0 = [
    ["B cell", "CL:0000236"],
    ["endothelial cell", "CL:0000115"],
    ["macrophage", "CL:0000235"],
    ["stromal cell", "CL:0000499"],
    ["erythroblast", "nan"],
    ["granulocyte", "nan"],
    ["basophil", "nan"],
    ["decidual stromal cell", "nan"],
    ["dendritic cell", "nan"],
    ["endodermal cell", "nan"],
    ["monocyte progenitor", "nan"],
    ["HSPC", "nan"],
    ["megakaryocte", "nan"],
    ["monocyte", "nan"],
    ["NK cell", "nan"],
    ["NKT cell", "nan"],
    ["PE lineage cell", "nan"],
    ["trophoblast progenitor", "nan"],
    ["labyrinthine trophoblast", "nan"],
    ["spiral artery trophoblast giant cells", "nan"],
    ["invasive spongiotrophoblast", "nan"],
    ["spongiotrophoblast", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_PLACENTA_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMousePlacenta(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_PLACENTA_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_PLACENTA_V0
        }
        super(CelltypeVersionsMousePlacenta, self).__init__(**kwargs)
