from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_MARROW_V0 = [
    ["basophil", "CL:0000767"],
    ["CD4-positive, alpha-beta T cell", "nan"],
    ["dendritic cell", "nan"],
    ["early pro-B cell", "CL:0002046"],
    ["erythroblast", "CL:0000765"],
    ["erythrocyte", "CL:0000232"],
    ["erythroid progenitor", "CL:0000038"],
    ["granulocyte monocyte progenitor cell", "nan"],
    ["granulocytopoietic cell", "CL:0002191"],
    ["hematopoietic precursor cell", "CL:0008001"],
    ["hematopoietic stem cell", "nan"],
    ["immature B cell", "CL:0000816"],
    ["late pro-B cell", "CL:0002048"],
    ["lymphoid progenitor cell", "nan"],
    ["macrophage", "nan"],
    ["mast cell", "nan"],
    ["monocyte", "CL:0000576"],
    ["megakaryocyte-erythroid progenitor cell", "CL:0000050"],
    ["naive B cell", "CL:0000788"],
    ["naive T cell", "CL:0000898"],
    ["neutrophil", "nan"],
    ["neutrophil progenitor", "nan"],
    ["NK cell", "CL:0000623"],
    ["plasma cell", "CL:0000786"],
    ["precursor B cell", "CL:0000817"],
    ["proerythroblast", "CL:0000547"],
    ["promonocyte", "CL:0000559"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_MARROW_V0 = {
    "names": {
        "granulocyte": ["basophil", "neutrophil", "mast cell"],
        "mature alpha-beta T cell": ["CD4-positive, alpha-beta T cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseMarrow(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_MARROW_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_MARROW_V0
        }
        super(CelltypeVersionsMouseMarrow, self).__init__(**kwargs)
