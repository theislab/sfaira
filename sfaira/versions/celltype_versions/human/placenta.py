from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_PLACENTA_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['CB CD34+', "nan"],
    ['Decidual Macrophages 1', "nan"],
    ['Decidual Macrophages 2', "nan"],
    ['Decidual Macrophages 3', "nan"],
    ['Decidual NK Cells 1', "nan"],
    ['Decidual NK Cells 2', "nan"],
    ['Decidual NK Cells 3', "nan"],
    ['Decidual NK Cells p', "nan"],
    ['Decidual Stromal Cells 1', "nan"],
    ['Decidual Stromal Cells 2', "nan"],
    ['Decidual Stromal Cells 3', "nan"],
    ['Dendritic Cells 1', "nan"],
    ['Dendritic Cells 2', "nan"],
    ['Endothelial Cells L', "nan"],
    ['Endothelial Cells f', "nan"],
    ['Endothelial Cells m', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epithelial Glandular Cells 1', "nan"],
    ['Epithelial Glandular Cells 2', "nan"],
    ['Epithelial cell (intermediated)', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Extravillous Trophoblasts', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal neuron', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fibroblasts 1', "nan"],
    ['Fibroblasts 2', "nan"],
    ['Granulocytes', "nan"],
    ['Hofbauer Cells', "nan"],
    ['ILC3', "nan"],
    ['Intermediated cell', "nan"],
    ['M2 Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Myeloid cell', "nan"],
    ['NK Cells CD16+', "nan"],
    ['NK Cells CD16-', "nan"],
    ['Neutrophil', "nan"],
    ['Perivascular Cells 1', "nan"],
    ['Perivascular Cells 2', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['Syncytiotrophoblasts', "nan"],
    ['T cell', "nan"],
    ['Villous Cytotrophoblasts', "nan"],
    ['hESC', "nan"]
]
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
            "0": CELLTYPES_HUMAN_PLACENTA_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_PLACENTA_V0
        }
        super(CelltypeVersionsHumanPlacenta, self).__init__(**kwargs)
