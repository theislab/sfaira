from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_PANCREAS_V0 = [
    ["b cell", "nan"],
    ["dendritic cell", "nan"],
    ["endothelial cell", "CL:0000115"],
    ["erythroblast", "nan"],
    ["glial cell", "nan"],
    ["granulocyte", "nan"],
    ["macrophage", "nan"],
    ["pancreatic A cell", "CL:0000171"],
    ["pancreatic acinar cell", "CL:0002064"],
    ["pancreatic B cell", "CL:0000169"],
    ["pancreatic D cell", "CL:0000173"],
    ["pancreatic ductal cell", "CL:0002079"],
    ["pancreatic PP cell", "CL:0002275"],
    ["pancreatic stellate cell", "CL:0002410"],
    ["smooth muscle cell", "nan"],
    ["stromal cell", "nan"],
    ["t cell", "nan"],
    ["lymphatic endothelial cell", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_PANCREAS_V0 = {
    "names": {
        "leukocyte": [
            "b cell",
            "dendritic cell",
            "granulocyte",
            "macrophage",
            "t cell"
        ],
    "endocrine cell": [
            "pancreatic A cell",
            "pancreatic D cell",
            "pancreatic PP cell"
        ],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMousePancreas(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_PANCREAS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_PANCREAS_V0
        }
        super(CelltypeVersionsMousePancreas, self).__init__(**kwargs)
