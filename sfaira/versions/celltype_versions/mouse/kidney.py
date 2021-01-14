from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_KIDNEY_V0 = [
    ["B cell", "CL:0000236"],
    ["brush cell", "nan"],
    ["dendritic cell", "nan"],
    ["endothelial cell", "nan"],
    ["fenestrated cell", "CL:0000666"],
    ["fetal adipocyte", "nan"],
    ["fetal mesenchymal cell", "nan"],
    ["fetal proliferative cell", "nan"],
    ["fibroblast", "CL:0000057"],
    ["interstitial fibroblast", "nan"],
    ["glomerular epithelial cell", "nan"],
    ["kidney collecting duct epithelial cell", "CL:1000454"],
    ["kidney collecting duct principal cell", "CL:1001431"],
    ["kidney cortex artery cell", "CL:1001045"],
    ["kidney distal convoluted tubule epithelial cell", "CL:1000849"],
    ["kidney loop of Henle ascending limb epithelial cell", "CL:1001016"],
    ["kidney loop of Henle thick ascending limb epithelial cell", "CL:1001106"],
    ["kidney proximal convoluted tubule epithelial cell", "CL:1000838"],
    ["kidney proximal straight tubule epithelial cell", "nan"],
    ["macrophage", "CL:0000235"],
    ["mesangial cell", "CL:0000650"],
    ["neutrophil progenitor", "nan"],
    ["NK cell", "nan"],
    ["podocyte", "CL:0000653"],
    ["plasma cell", "CL:0000786"],
    ["T cell", "CL:0000084"],
    ["ureteric epithelial cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_KIDNEY_V0 = {
    "names": {
        "epithelial cell": [
            "kidney collecting duct epithelial cell",
            "kidney distal convoluted tubule epithelial cell",
            "kidney loop of Henle ascending limb epithelial cell",
            "kidney loop of Henle thick ascending limb epithelial cell",
            "kidney proximal convoluted tubule epithelial cell",
            "kidney proximal straight tubule epithelial cell",
        ],
        "epithelial cell of proximal tubule": [
            "kidney proximal convoluted tubule epithelial cell",
            "kidney proximal straight tubule epithelial cell",
        ],
        "lymphocyte": ["B cell", "dendritic cell", "macrophage", "NK cell", "T cell"],
        "leukocyte": [
            "B cell", "dendritic cell", "macrophage", "neutrophil progenitor",
            "NK cell", "plasma cell", "T cell"
        ],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseKidney(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_KIDNEY_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_KIDNEY_V0
        }
        super(CelltypeVersionsMouseKidney, self).__init__(**kwargs)
