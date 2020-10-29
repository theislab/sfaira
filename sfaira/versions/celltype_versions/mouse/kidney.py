from .external import CelltypeVersionsBase

# Version 0
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
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_KIDNEY_V0
        }
        super(CelltypeVersionsMouseKidney, self).__init__(**kwargs)
