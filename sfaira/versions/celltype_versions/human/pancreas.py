from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_PANCREAS_V0 = [
    ['Acinar cell', "nan"],
    ['Activated Stellate cell', "nan"],
    ['Alpha cell', "nan"],
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['Beta cell', "nan"],
    ['CB CD34+', "nan"],
    ['Co-expression cell', "nan"],
    ['Delta cell', "nan"],
    ['Dendritic cell', "nan"],
    ['Ductal cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Enterocyte', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Epithelial progenitor', "nan"],
    ['Epsilon cell', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fibroblast', "nan"],
    ['Gamma cell', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Immature sertoli cell (Pre-Sertoli cell)', "nan"],
    ['MHC class II cell', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Mesenchymal Cell', "nan"],
    ['Monocyte', "nan"],
    ['Neuron', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['PSC cell', "nan"],
    ['Pancreas exocrine cell', "nan"],
    ['Primordial germ cell', "nan"],
    ['Proximal tubule progenitor', "nan"],
    ['Quiescent Stellate cell', "nan"],
    ['Schwann cell', "nan"],
    ['Skeletal muscle cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['Unclassified endocrine cell', "nan"],
    ['Unknown', "nan"]
]
ONTOLOGIES_HUMAN_PANCREAS_V0 = {
    "names": {
        'Endocrine cell': ['Alpha cell', 'Beta cell', 'Gamma cell', 'Delta cell', 'Epsilon cell', 'Unclassified endocrine cell'],
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanPancreas(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_PANCREAS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_PANCREAS_V0
        }
        super(CelltypeVersionsHumanPancreas, self).__init__(**kwargs)
