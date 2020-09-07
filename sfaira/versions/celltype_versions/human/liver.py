from .external import CelltypeVersionsBase

CELLTYPES_HUMAN_LIVER_V0 = [
    ['Alpha beta T cells', "nan"],
    ['Antigen presenting cell (RPS high)', "nan"],
    ['CB CD34+', "nan"],
    ['Central venous LSECs', "nan"],
    ['Cholangiocytes', "nan"],
    ['Dendritic cell 1', "nan"],
    ['Dendritic cell 2', "nan"],
    ['Dendritic cell precursor', "nan"],
    ['Early Erythroid', "nan"],
    ['Early lymphoid T lymphocyte', "nan"],
    ['Endothelial cell (APC)', "nan"],
    ['Endothelial cell (endothelial to mesenchymal transition)', "nan"],
    ['Enterocyte ', "nan"],
    ['Enterocyte progenitor', "nan"],
    ['Epithelial progenitor', "nan"],
    ['Fibroblast', "nan"],
    ['Gamma delta T cells 1', "nan"],
    ['Gamma delta T cells 2', "nan"],
    ['Gastric endocrine cell', "nan"],
    ['Goblet cell', "nan"],
    ['HSC MPP', "nan"],
    ['Hepatic stellate cells', "nan"],
    ['Hepatocyte 1', "nan"],
    ['Hepatocyte 2', "nan"],
    ['Hepatocyte 3', "nan"],
    ['Hepatocyte 4', "nan"],
    ['Hepatocyte 5', "nan"],
    ['Hepatocyte 6', "nan"],
    ['ILC', "nan"],
    ['ILC precursor', "nan"],
    ['Inflammatory macrophages', "nan"],
    ['Kupffer Cell', "nan"],
    ['Late Erythroid', "nan"],
    ['Liver sinusoidal endothelial cells', "nan"],
    ['MEMP', "nan"],
    ['MP', "nan"],
    ['Macrovascular endothelial cells', "nan"],
    ['Mast cell', "nan"],
    ['Mature B cells', "nan"],
    ['Megakaryocyte', "nan"],
    ['Mesenchyme', "nan"],
    ['Mesothelia', "nan"],
    ['Mid Erythroid', "nan"],
    ['Mono Macrophage', "nan"],
    ['Monocyte', "nan"],
    ['Monocyte precursor', "nan"],
    ['Myeloid cell', "nan"],
    ['NK cell', "nan"],
    ['Neutrophil', "nan"],
    ['Neutrophil (RPS high)', "nan"],
    ['Neutrophil myeloid progenitor', "nan"],
    ['Non inflammatory macrophages', "nan"],
    ['Other endothelial cells', "nan"],
    ['Pancreas exocrine cell', "nan"],
    ['Periportal LSECs', "nan"],
    ['Plasma B cell', "nan"],
    ['Plasma cells', "nan"],
    ['Pre pro B cell', "nan"],
    ['Primordial germ cell', "nan"],
    ['Proliferating T cell', "nan"],
    ['Smooth muscle cell', "nan"],
    ['Unknown', "nan"],
    ['VCAM1pos EI macrophage', "nan"],
    ['pDendritic cell precursor', "nan"],
    ['pre B cell', "nan"],
    ['pro B cell', "nan"]
]
ONTOLOGIES_HUMAN_LIVER_V0 = {
    "names": {
        'Erythroid cells': ['Early Erythroid', 'Mid Erythroid', 'Late Erythroid'],
        'Endothelial cell': ['Liver sinusoidal endothelial cells', 'Macrovascular endothelial cells', 'Other endothelial cells'],
        'Hepatocyte': ['Hepatocyte 1','Hepatocyte 2','Hepatocyte 3','Hepatocyte 4','Hepatocyte 5','Hepatocyte 6'],
        'Hepatocytes': ['Hepatocyte 1','Hepatocyte 2','Hepatocyte 3','Hepatocyte 4','Hepatocyte 5','Hepatocyte 6'],
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
            "0": CELLTYPES_HUMAN_LIVER_V0
        }
        self.ontology = {
            "0": ONTOLOGIES_HUMAN_LIVER_V0
        }
        super(CelltypeVersionsHumanLiver, self).__init__(**kwargs)
