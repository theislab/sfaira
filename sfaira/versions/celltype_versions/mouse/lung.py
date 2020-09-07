from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_LUNG_V0 = [
    ["adventitial cell", "nan"],
    ["alveolar bipotent progenitor", "nan"],
    ["alveolar epithelial cell type I", "nan"],
    ["alveolar epithelial cell type II", "nan"],
    ["alveolar macrophage", "nan"],
    ["B cell", "CL_0000236"],
    ["basophil", "nan"],
    ["bronchial smooth muscle cell", "nan"],
    ["CD4-positive, alpha-beta T cell", "nan"],
    ["CD8-positive, alpha-beta T cell", "nan"],
    ["ciliated cell", "nan"],
    ["clara cell", "nan"],
    ["classical monocyte", "nan"],
    ["club cell of bronchiole", "nan"],
    ["endothelial cell of lymphatic vessel", "nan"],
    ["eosinophil", "nan"],
    ["fibroblast of lung", "nan"],
    ["glial cell", "CL_0000125"],
    ["intermediate monocyte", "nan"],
    ["lung macrophage", "nan"],
    ["lung neuroendocrine cell", "nan"],
    ["monocyte progenitor", "nan"],
    ["myeloid dendritic cell", "nan"],
    ["neutrophil", "nan"],
    ["NK cell", "CL_0000623"],
    ["non-classical monocyte", "nan"],
    ["nuocyte", "nan"],
    ["pericyte cell", "nan"],
    ["plasma cell", "nan"],
    ["plasmacytoid dendritic cell", "nan"],
    ["proliferative cell", "nan"],
    ["pulmonary interstitial fibroblast", "nan"],
    ["regulatory T cell", "nan"],
    ["respiratory basal cell", "nan"],
    ["smooth muscle cell of the pulmonary artery", "nan"],
    ["type I pneumocyte", "nan"],
    ["type II pneumocyte", "nan"],
    ["vein endothelial cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_LUNG_V0 = {
    "names": {
        "dendritic cell": [
            "myeloid dendritic cell", "plasmacytoid dendritic cell"
        ],
        "endothelial cell": [
            "endothelial cell of lymphatic vessel", "vein endothelial cell"
        ],
        "leukocyte": [
            "alveolar macrophage", "B cell", "basophil", "bronchial smooth muscle cell",
            "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "classical monocyte",
            "eosinophil", "glial cell", "intermediate monocyte", "lung macrophage", "monocyte progenitor",
            "myeloid dendritic cell", "neutrophil", "NK cell", "non-classical monocyte", "plasma cell",
            "plasmacytoid dendritic cell", "regulatory T cell"
        ],
        "lymphocyte": [
            "B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
            "plasmacytoid dendritic cell", "NK cell", "regulatory T cell"
        ],
        "mature NK T cell": [
            "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "NK cell", "regulatory T cell"
        ],
        "T cell": [
            "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "regulatory T cell"
        ],
        "stromal cell": [
            "fibroblast of lung", "pericyte cell"
        ]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseLung(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_LUNG_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LUNG_V0
        }
        super(CelltypeVersionsMouseLung, self).__init__(**kwargs)
