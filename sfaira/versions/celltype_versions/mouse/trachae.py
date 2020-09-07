from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_TRACHAE_V0 = [
    ["basal epithelial cell of tracheobronchial tree", "CL:0002329"],
    ["chondrocyte", "CL:0000138"],
    ["ciliated columnar cell of tracheobronchial tree", "CL:0002145"],
    ["endothelial cell", "CL:000115"],
    ["epithelial cell", "CL:000115"],
    ["fibroblast", "CL:0000057"],
    ["granulocyte", "CL:0000094"],
    ["keratinocyte", "nan"],
    ["macrophage", "CL:0000235"],
    ["mesenchymal cell", "nan"],
    ["mesenchymal progenitor cell", "nan"],
    ["mucus secreting cell", "CL:0000319"],
    ["neuroendocrine cell", "nan"],
    ["smooth muscle cell of trachea", "CL:0002600"],
    ["T cell", "CL:0000084"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_TRACHAE_V0 = {
    "names": {
        'blood cell': ["granulocyte", "macrophage", "T cell"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseTrachae(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_TRACHAE_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_TRACHAE_V0
        }
        super(CelltypeVersionsMouseTrachae, self).__init__(**kwargs)
