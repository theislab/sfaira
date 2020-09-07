from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_LARGEINTESTINE_V0 = [
    ["Brush cell of epithelium proper of large intestine", "CL:0002203"],
    ["enterocyte of epithelium of large intestine", "CL:0002071"],
    ["enteroendocrine cell", "CL:0000164"],
    ["epithelial cell of large intestine", "CL:0002253"],
    ["hematopoietic stem cell", "CL:0000037"],
    ["intestinal crypt stem cell", "CL:0002250"],
    ["large intestine goblet cell", "CL:1000320"],
    ["secretory cell", "CL:0000151"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_LARGEINTESTINE_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseLargeintestine(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_LARGEINTESTINE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LARGEINTESTINE_V0
        }
        super(CelltypeVersionsMouseLargeintestine, self).__init__(**kwargs)
