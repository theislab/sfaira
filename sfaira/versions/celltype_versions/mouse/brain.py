from .external import CelltypeVersionsBase

# Version 0
CELLTYPES_MOUSE_BRAIN_V0 = [
    ["astrocyte", "CL:0000127"],
    ["BAM", "nan"],
    ["B cells", "nan"],
    ["Bergmann glial cell", "CL:0000644"],
    ["brain pericyte", "CL:2000043"],
    ["CD8-positive, alpha-beta T cell", "CL:0000625"],
    ["cDC1", "nan"],
    ["cDC2", "nan"],
    ["endothelial cell", "CL:0000115"],
    ["ependymal cell", "CL:0000065"],
    ["GABAergic cell", "nan"],
    ["granulocyte", "nan"],
    ["ILC", "nan"],
    ["interneuron", "CL:0000099"],
    ["macrophage", "CL:0000235"],
    ["mature NK T cell", "nan"],
    ["medium spiny neuron", "CL:1001474"],
    ["microglial cell", "CL:0000129"],
    ["migDC", "nan"],
    ["monocyte", "nan"],
    ["neuroepithelial cell", "nan"],
    ["neuron", "CL:0000540"],
    ["neuronal stem cell", "CL:0000047"],
    ["neutorphils", "nan"],
    ["NK cells", "nan"],
    ["oligodendrocyte", "CL:0000128"],
    ["oligodendrocyte precursor cell", "CL:0002453"],
    ["pDC", "nan"],
    ["schwann cell", "nan"],
    ["yd T cells", "nan"],
    ["unknown", "nan"]
]
ONTOLOGIES_MOUSE_BRAIN_V0 = {
    "names": {
        "T cell": ["CD8-positive, alpha-beta T cell", "yd T cells", "mature NK T cell"],
    },
    "ontology_ids": {},
}


class CelltypeVersionsMouseBrain(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_MOUSE_BRAIN_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_MOUSE_BRAIN_V0
        }
        super(CelltypeVersionsMouseBrain, self).__init__(**kwargs)
