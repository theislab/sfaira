from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_ESOPHAGUS_V0 = [
    ['B cell (Plasmocyte)', "nan"],
    ['B_CD27neg', "nan"],
    ['B_CD27pos', "nan"],
    ['Basal cell', "nan"],
    ['Blood_vessel', "nan"],
    ['CB CD34+', "nan"],
    ['Dendritic cell', "nan"],
    ['Endothelial cell', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Epi_dividing', "nan"],
    ['Epi_suprabasal', "nan"],
    ['Epi_upper', "nan"],
    ['Erythroid progenitor cell (RP high)', "nan"],
    ['Fetal epithelial progenitor', "nan"],
    ['Fetal mesenchymal progenitor', "nan"],
    ['Fetal stromal cell', "nan"],
    ['Fibroblast', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Glands_duct', "nan"],
    ['Glands_mucous', "nan"],
    ['Loop of Henle', "nan"],
    ['Lymph_vessel', "nan"],
    ['Macrophage', "nan"],
    ['Mast cell', "nan"],
    ['Monocyte', "nan"],
    ['NK_T_CD8_Cytotoxic', "nan"],
    ['Neutrophil', "nan"],
    ['Sinusoidal endothelial cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Stratified epithelial cell', "nan"],
    ['Stromal cell', "nan"],
    ['T_CD4', "nan"],
    ['T_CD8', "nan"]
]
ONTOLOGIES_HUMAN_ESOPHAGUS_V0 = {
    "names": {
        "Mono_macro": ["Monocyte", "Macrophage"],
        "B cell": ['B_CD27neg', 'B_CD27pos'],
        "T cell": ["T_CD4", "T_CD8", "NK_T_CD8_Cytotoxic"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanEsophagus(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_ESOPHAGUS_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_ESOPHAGUS_V0
        }
        super(CelltypeVersionsHumanEsophagus, self).__init__(**kwargs)
