from .external import CelltypeVersionsBase

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
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_COLON_V0
        }
        super(CelltypeVersionsHumanColon, self).__init__(**kwargs)
