from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_EYE_V0 = [
    ['Amacrine cell', "nan"],
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B-cell', "nan"],
    ['Basal cell', "nan"],
    ['CB CD34_pos', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Epithelial cell (intermediated)', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal endocrine cell', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal neuron', "nan"],
    ['Fetal skeletal muscle cell', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Goblet cell', "nan"],
    ['Horizontal cells', "nan"],
    ['Macroglia', "nan"],
    ['Macrophage', "nan"],
    ['Mast-cell', "nan"],
    ['Melanocyte', "nan"],
    ['Microglia', "nan"],
    ['Muller cell', "nan"],
    ['Pericyte', "nan"],
    ['Primordial germ cell', "nan"],
    ['Retinal bipolar neuron type A', "nan"],
    ['Retinal bipolar neuron type B', "nan"],
    ['Retinal bipolar neuron type C', "nan"],
    ['Retinal bipolar neuron type D', "nan"],
    ['Retinal cone cell', "nan"],
    ['Retinal ganglion cell', "nan"],
    ['Retinal pigment epithelium', "nan"],
    ['Retinal rod cell type A', "nan"],
    ['Retinal rod cell type B', "nan"],
    ['Retinal rod cell type C', "nan"],
    ['Schwann1', "nan"],
    ['Schwann2', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['T cell', "nan"],
    ['T/NK-cell', "nan"],
    ['Unknown', "nan"],
    ['hESC', "nan"]
]
ONTOLOGIES_HUMAN_EYE_V0 = {
    "names": {
        'BPs': ['Retinal bipolar neuron type A', 'Retinal bipolar neuron type B', 'Retinal bipolar neuron type C', 'Retinal bipolar neuron type D'],
        'Rods': ['Retinal rod cell type A', 'Retinal rod cell type B', 'Retinal rod cell type C', ]
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanEye(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_EYE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_EYE_V0
        }
        super(CelltypeVersionsHumanEye, self).__init__(**kwargs)
