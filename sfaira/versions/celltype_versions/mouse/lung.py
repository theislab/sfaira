from .external import CelltypeVersionsBase

# Version 0
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
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_LUNG_V0
        }
        super(CelltypeVersionsMouseLung, self).__init__(**kwargs)
