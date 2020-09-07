from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_HEART_V0 = [
    ["aortic endothelial cell", "CL:0002544"],
    ["atrial myocyte", "CL:0002129"],
    ["B cell", "CL:CL:0000115"],
    ["cardiac neuron", "CL:0000057"],
    ["cardiomyocyte", "CL:0000746"],
    ["endocardial cell", "CL:0002350"],
    ["endothelial cell of coronary artery", "CL:2000018"],
    ["epithelial cell", "CL:"],
    ["erythrocyte", "CL:"],
    ["fibroblast of cardiac tissue", "CL:0002548"],
    ["fibrocyte", "CL:CL:0000145"],
    ["leukocyte", "CL:0000738"],
    ["mast cell", "nan"],
    ["monocyte", "nan"],
    ["macrophage", "nan"],
    ["professional antigen presenting cell", "nan"],
    ["smooth muscle cell", "CL:0000192"],
    ["T cell", "nan"],
    ["valve cell", "CL:0000663"],
    ["ventricular myocyte", "CL:0002131"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_HEART_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsMouseHeart(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_HEART_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_HEART_V0
        }
        super(CelltypeVersionsMouseHeart, self).__init__(**kwargs)
