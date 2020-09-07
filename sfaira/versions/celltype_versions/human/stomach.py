from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_STOMACH_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['Basal cell', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Enterocyte', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Epithelial cell', "nan"],
    ['Epithelial cell (intermediated)', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fasciculata cell', "nan"],
    ['Fetal Neuron', "nan"],
    ['Fetal acinar cell', "nan"],
    ['Fetal chondrocyte', "nan"],
    ['Fetal endocrine cell', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal fibroblast', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal neuron', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric chief cell', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Goblet cell', "nan"],
    ['Hepatocyte/Endodermal cell', "nan"],
    ['M2 Macrophage', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Mesothelial cell', "nan"],
    ['Monocyte', "nan"],
    ['Myeloid cell', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Primordial germ cell', "nan"],
    ['Proliferating T cell', "nan"],
    ['Proximal tubule progenitor', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['hESC', "nan"]
]
ONTOLOGIES_HUMAN_STOMACH_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanStomach(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_STOMACH_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_STOMACH_V0
        }
        super(CelltypeVersionsHumanStomach, self).__init__(**kwargs)
