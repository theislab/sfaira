from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_HESC_V0 = [
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal neuron', "nan"],
    ['Primordial germ cell', "nan"],
    ['Proliferating T cell', "nan"],
    ['hESC', "nan"]
]
ONTOLOGIES_HUMAN_HESC_V0 = {
    "names": {},
    "ontology_ids": {},
}


class CelltypeVersionsHumanHesc(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_HESC_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_HESC_V0
        }
        super(CelltypeVersionsHumanHesc, self).__init__(**kwargs)
