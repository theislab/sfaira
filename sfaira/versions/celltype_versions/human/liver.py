from .external import CelltypeVersionsBase

ONTOLOGIES_HUMAN_LIVER_V0 = {
    "names": {
        'Erythroid cells': ['Early Erythroid', 'Mid Erythroid', 'Late Erythroid'],
        'Endothelial cell': ['Liver sinusoidal endothelial cells', 'Macrovascular endothelial cells', 'Other endothelial cells'],
        'Hepatocyte': ['Hepatocyte 1', 'Hepatocyte 2', 'Hepatocyte 3', 'Hepatocyte 4', 'Hepatocyte 5', 'Hepatocyte 6'],
        'Hepatocytes': ['Hepatocyte 1', 'Hepatocyte 2', 'Hepatocyte 3', 'Hepatocyte 4', 'Hepatocyte 5', 'Hepatocyte 6'],
        'Endothelia': ['Liver sinusoidal endothelial cells', 'Macrovascular endothelial cells', 'Other endothelial cells'],
        'Bcells': ['pro B cell', 'Pre pro B cell', 'Mature B cells', 'pre B cell', 'Plasma B cell'],
        'Tcells': ['Gamma delta T cells 2', 'Gamma delta T cells 1', 'Alpha beta T cells'],
        'pDCs': ['Dendritic cell 1', 'Dendritic cell 2'],
        'NK, NKT and T cells': ['NK cell', 'Alpha beta T cells', 'Gamma delta T cells 1', 'Gamma delta T cells 2'],
        'B Cell': ['pro B cell', 'Pre pro B cell', 'Mature B cells', 'pre B cell', 'Plasma B cell'],
        'T cell': ['Alpha beta T cells', 'Gamma delta T cells 1', 'Gamma delta T cells 2'],
        'Dendritic cell': ['Dendritic cell 1', 'Dendritic cell 2'],
        'B cell': ['pro B cell', 'Pre pro B cell', 'Mature B cells', 'pre B cell']
    },
    "ontology_ids": {},
}


class CelltypeVersionsHumanLiver(CelltypeVersionsBase):

    def __init__(self, **kwargs):
        self.celltype_universe = {
            "0": self.read_csv(".".join(__file__.split(".")[:-1])+".csv")
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_LIVER_V0
        }
        super(CelltypeVersionsHumanLiver, self).__init__(**kwargs)
