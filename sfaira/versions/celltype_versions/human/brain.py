from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_BRAIN_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['Astrocytes 1', "nan"],
    ['Astrocytes 2', "nan"],
    ['B cell', "nan"],
    ['B cell (Plasmocyte)', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cells', "nan"],
    ['Erythroid cell', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal endocrine cell', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['GABAergic interneurons 1', "nan"],
    ['GABAergic interneurons 2', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Glutamatergic neurons from the PFC 1', "nan"],
    ['Glutamatergic neurons from the PFC 2', "nan"],
    ['Goblet cell', "nan"],
    ['Granule neurons from the hip dentate gyrus region', "nan"],
    ['Macrophage', "nan"],
    ['Microglia', "nan"],
    ['Monocyte', "nan"],
    ['Neuronal stem cells', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Oligodendrocyte precursors', "nan"],
    ['Oligodendrocytes', "nan"],
    ['Primordial germ cell', "nan"],
    ['Pyramidal neurons from the hip CA region 1', "nan"],
    ['Pyramidal neurons from the hip CA region 2', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stromal cell', "nan"],
    ['T cell', "nan"],
    ['Unknown', "nan"]
]
ONTOLOGIES_HUMAN_BRAIN_V0 = {
    "names": {
        'Astrocyte': ['Astrocytes 1', 'Astrocytes 2'],
        'Fetal Neuron': ['Glutamatergic neurons from the PFC 1', 'Glutamatergic neurons from the PFC 2',
                         'Granule neurons from the hip dentate gyrus region', 'GABAergic interneurons 1',
                         'GABAergic interneurons 2', 'Pyramidal neurons from the hip CA region 1', 'Pyramidal neurons from the hip CA region 2']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanBrain(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_BRAIN_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_BRAIN_V0
        }
        super(CelltypeVersionsHumanBrain, self).__init__(**kwargs)
