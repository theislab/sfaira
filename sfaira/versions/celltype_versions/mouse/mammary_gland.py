from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_MAMMARYGLAND_V0 = [
    ["B cell", "CL:0000236"],
    ["basal cell", "CL:0000646"],
    ["endothelial cell", "CL:0000115"],
    ["luminal epithelial cell of mammary gland", "CL:0002326"],
    ["luminal progenitor cell", "CL:0002326"],
    ["macrophage", "CL:0000235"],
    ["stromal cell", "CL:0000499"],
    ["T cell", "CL:0000084"],
    ["dendritic cell", "nan"],
    ["proliferative cell", "nan"],
    ["NK cell", "CL:0000623"],
    ["stem and progenitor cell", "CL:0000623"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_MAMMARYGLAND_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseMammarygland(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_MAMMARYGLAND_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_MAMMARYGLAND_V0
        }
        super(CelltypeVersionsMouseMammarygland, self).__init__(**kwargs)
