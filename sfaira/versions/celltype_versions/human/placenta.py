from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_PLACENTA_V0 = {
    "names": {
        'Fibroblast': ['Fibroblasts 1', 'Fibroblasts 2'],
        'Macrophage': ['Decidual Macrophages 1', 'Decidual Macrophages 2', 'Decidual Macrophages 3'],
        'Epithelial cell': ['Epithelial Glandular Cells 1', 'Epithelial Glandular Cells 2'],
        'Fetal stromal cell': ['Decidual Stromal Cells 1', 'Decidual Stromal Cells 2', 'Decidual Stromal Cells 3'],
        'Endothelial cell': ['Endothelial Cells f', 'Endothelial Cells m', 'Endothelial Cells L'],
        'Dendritic cell': ['Dendritic Cells 1', 'Dendritic Cells 2'],
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanPlacenta(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_PLACENTA_V0
        }
        super(CelltypeVersionsHumanPlacenta, self).__init__(**kwargs)
