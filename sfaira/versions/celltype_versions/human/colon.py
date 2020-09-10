from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_COLON_V0 = [
    ['Antigen presenting cell (RPS high)', "nan"],
    ['B cell IgA Plasma', "nan"],
    ['B cell IgG Plasma', "nan"],
    ['B cell cycling', "nan"],
    ['B cell memory', "nan"],
    ['Best4+ Enterocytes', "nan"],
    ['CD4+ Memory', "nan"],
    ['CD4+ PD1+', "nan"],
    ['CD4+ T Activated Fos-hi', "nan"],
    ['CD4+ T Activated Fos-lo', "nan"],
    ['CD69+ Mast', "nan"],
    ['CD69- Mast', "nan"],
    ['CD8 T', "nan"],
    ['CD8+ IELs', "nan"],
    ['CD8+ IL17+', "nan"],
    ['CD8+ LP', "nan"],
    ['Cycling T', "nan"],
    ['Cycling TA', "nan"],
    ['DC1', "nan"],
    ['DC2', "nan"],
    ['Endothelial', "nan"],
    ['Enterocyte Progenitors', "nan"],
    ['Enterocytes', "nan"],
    ['Enteroendocrine cells', "nan"],
    ['Erythroid cell', "nan"],
    ['Fetal Neuron', "nan"],
    ['Fetal enterocyte ', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fibroblast', "nan"],
    ['Follicular', "nan"],
    ['Glial cells', "nan"],
    ['Goblet cells', "nan"],
    ['ILC', "nan"],
    ['Immature Enterocytes 1', "nan"],
    ['Immature Enterocytes 2', "nan"],
    ['Immature Goblet', "nan"],
    ['LYVE1 Macrophage', "nan"],
    ['Lymphoid DC', "nan"],
    ['M cells', "nan"],
    ['MT-hi', "nan"],
    ['Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Myofibroblasts', "nan"],
    ['NK', "nan"],
    ['Neutrophil', "nan"],
    ['Paneth cells', "nan"],
    ['Pericytes', "nan"],
    ['Primordial germ cell', "nan"],
    ['Secretory TA', "nan"],
    ['Smooth Muscle', "nan"],
    ['Stem cells', "nan"],
    ['Stromal', "nan"],
    ['TA 1', "nan"],
    ['TA 2', "nan"],
    ['Tcm', "nan"],
    ['Tfh', "nan"],
    ['Th1', "nan"],
    ['Th17', "nan"],
    ['Treg', "nan"],
    ['Tregs', "nan"],
    ['Tuft', "nan"],
    ['WNT2B+ Fos-lo 1', "nan"],
    ['WNT5B+ 2', "nan"],
    ['cycling DCs', "nan"],
    ['cycling gd T', "nan"],
    ['gd T', "nan"],
    ['pDC', "nan"]
]
ONTOLOGIES_HUMAN_COLON_V0 = {
    "names": {
        'Plasma Cells': ['B cell IgA Plasma', 'B cell IgG Plasma'],
        'Macrophage': ['LYVE1 Macrophage', 'Macrophage'],
        'Enterocytes': ['Enterocytes', 'Best4+ Enterocytes'],
        'TA': ['Cycling TA', 'TA 1', 'TA 2', 'Secretory TA'],
        'Activated CD4 T': ['CD4+ T Activated Fos-hi', 'CD4+ T Activated Fos-lo'],
        'Fetal enterocyte': ['Immature Enterocytes 1', 'Immature Enterocytes 2'],
        'B cell (Plasmocyte)': ['B cell IgA Plasma', 'B cell IgG Plasma'],
        'Mast cell': ['CD69+ Mast', 'CD69- Mast'],
        'Dendritic cell': ['DC1', 'DC2'],
        'B cell': ['B cell cycling', 'B cell memory', 'Follicular'],
        'T cell': ['Treg', 'Cycling T', 'CD4+ T Activated Fos-hi', 'CD4+ T Activated Fos-lo', 'Tcm', 'Tfh', 'Th1', 'Th17', 'cycling gd T', 'gd T'],
        'Epithelial cell': ['Enterocytes', 'Goblet cells', 'Enteroendocrine cells', 'Tuft']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanColon(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_COLON_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_COLON_V0
        }
        super(CelltypeVersionsHumanColon, self).__init__(**kwargs)
