from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_KIDNEY_V0 = {
    "names": {
        'Type A intercalated cell': ['Collecting Duct - Intercalated Cells Type A (cortex)',
                                         'Collecting Duct - Intercalated Cells Type A (medulla)'],
        'Principal cell': ['Collecting Duct - PCs - Stressed Dissoc Subset',
                           'Collecting Duct - Principal Cells (cortex)',
                           'Collecting Duct - Principal Cells (medulla)'],
        'Proximal tubule': ['Proximal Tubule Epithelial Cells (S1)',
                            'Proximal Tubule Epithelial Cells (S2)',
                            'Proximal Tubule Epithelial Cells (S3)',
                            'Proximal Tubule Epithelial Cells - Fibrinogen+ (S3)',
                            'Proximal Tubule Epithelial Cells - Stress/Inflam'],
        'Dendritic cell': ['MNP-c/dendritic cell', 'Plasmacytoid dendritic cell'],
        'Endothelial cell': ['Endothelial Cells (unassigned)', 'Endothelial Cells - AEA & DVR', 'Endothelial Cells - AVR', 'Endothelial Cells - glomerular capillaries', 'Peritubular capillary endothelium 1', 'Peritubular capillary endothelium 2'],
        'Epithelial cell': ['Pelvic epithelium', 'Pelvic epithelium - distal UB', 'Proximal Tubule Epithelial Cells (S1)', 'Proximal Tubule Epithelial Cells (S2)', 'Proximal Tubule Epithelial Cells (S3)', 'Proximal Tubule Epithelial Cells - Fibrinogen+ (S3)', 'Proximal Tubule Epithelial Cells - Stress/Inflam'],
        'Intercalated cell': ['Collecting Duct - Intercalated Cells Type A (cortex)', 'Collecting Duct - Intercalated Cells Type A (medulla)', 'Collecting Duct - Intercalated Cells Type B', 'Indistinct intercalated cell'],
        'T cell': ['CD4 T cell', 'CD8 T cell'],
        'Ureteric bud cell': ['CNT/PC - proximal UB', 'Proximal UB', 'Pelvic epithelium - distal UB']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanKidney(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_KIDNEY_V0
        }
        super(CelltypeVersionsHumanKidney, self).__init__(**kwargs)
