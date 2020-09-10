from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_LUNG_V0 = [
    ['Cycling cells', "nan"],
    ['Cartilage', "nan"],
    ['Fetal airway progenitors', "nan"],
    ['Mesothelium', "nan"],
    ['AT1', "nan"],
    ['AT2', "nan"],
    ['Acinar', "nan"],
    ['Airway smooth muscle', "nan"],
    ['Arterial', "nan"],
    ['B cell lineage', "nan"],
    ['Basal', "nan"],
    ['Bronchial Vessel 1', "nan"],
    ['Bronchial Vessel 2', "nan"],
    ['Capillary', "nan"],
    ['Capillary Intermediate 1', "nan"],
    ['Capillary Intermediate 2', "nan"],
    ['Dendritic cells', "nan"],
    ['Erythrocytes', "nan"],
    ['Fibroblasts', "nan"],
    ['Fibromyocyte', "nan"],
    ['Innate lymphoid cells', "nan"],
    ['KRT5-/KRT17+', "nan"],
    ['Lymphatic EC', "nan"],
    ['Macrophages', "nan"],
    ['Mast cells', "nan"],
    ['Megakaryocytes', "nan"],
    ['Monocytes', "nan"],
    ['Multiciliated lineage', "nan"],
    ['Myofibroblasts', "nan"],
    ['Neutrophilic', "nan"],
    ['Proliferating Epithelial Cells', "nan"],
    ['Rare', "nan"],
    ['Secretory', "nan"],
    ['Submucosal Secretory', "nan"],
    ['T cell lineage', "nan"],
    ['Venous', "nan"],
    ['Venous smooth muscle', "nan"],
    ['unknown', "nan"]
]
ONTOLOGIES_HUMAN_LUNG_V0 = {
    "names": {
        "1_Endothelial": ['Arterial', 'Capillary', 'Venous', 'Bronchial Vessel 1', 'Bronchial Vessel 2',
                          'Capillary Intermediate 1', 'Capillary Intermediate 2', 'Lymphatic EC'],
        "1_Epithelial": ['Basal', 'Multiciliated lineage', 'Secretory', 'Rare', 'Submucosal Secretory', 'Acinar',
                         'AT1', 'AT2', 'KRT5-/KRT17+', 'Proliferating Epithelial Cells',
                         'Fetal airway progenitors'],
        "1_Immune": ["B cell lineage", "T cell lineage", "Innate lymphoid cells", "Dendritic cells", "Macrophages",
                     "Monocytes", "Mast cells", "Megakaryocytes", "Erythrocytes"],
        "1_Stroma": ['Mesothelium', 'Fibroblasts', 'Myofibroblasts', 'Fibromyocyte', 'Airway smooth muscle',
                     'Venous smooth muscle', 'Cartilage'],
        "2_Blood vessels": ['Arterial', 'Capillary', 'Venous', 'Bronchial Vessel 1', 'Bronchial Vessel 2',
                            'Capillary Intermediate 1', 'Capillary Intermediate 2'],
        "2_Fibroblast lineage": ['Fibroblasts', 'Myofibroblasts'],
        "2_Lymphoid": ['B cell lineage', 'T cell lineage', 'Innate lymphoid cells'],
        "2_Smooth Muscle": ['Fibromyocyte', 'Airway smooth muscle', 'Venous smooth muscle'],
        "2_Myeloid": ["Dendritic cells", "Macrophages", "Monocytes", "Mast cells"]
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanLung(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": CELLTYPES_HUMAN_LUNG_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_LUNG_V0
        }
        super(CelltypeVersionsHumanLung, self).__init__(**kwargs)
