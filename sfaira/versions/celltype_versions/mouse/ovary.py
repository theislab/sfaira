from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_OVARY_V0 = [
    ["cumulus cell", "nan"],
    ["granulosa cell", "nan"],
    ["large luteal cell", "nan"],
    ["macrophage", "nan"],
    ["small luteal cell", "nan"],
    ["epithelial cell of ovarian surface", "nan"],
    ["endothelial cell of ovarian surface", "nan"],
    ["stromal cell", "nan"],
    ["thecal cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_OVARY_V0 = {
    "names": {
        'luteal cell': ['small luteal cell', 'large luteal cell'],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseOvary(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_OVARY_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_OVARY_V0
        }
        super(CelltypeVersionsMouseOvary, self).__init__(**kwargs)
