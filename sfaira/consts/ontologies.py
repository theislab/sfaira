from sfaira.versions.metadata import OntologyList, OntologyCelltypes
from sfaira.versions.metadata import OntologyUberon, OntologyHsapdv, OntologyMmusdv, \
    OntologySinglecellLibraryConstruction, OntologyCellosaurus


class OntologyContainerSfaira:

    _cellontology_class: OntologyCelltypes

    def __init__(self):
        self.age = None
        self.annotated = OntologyList(terms=[True, False])
        self.author = None
        self.assay_differentiation = None
        self.assay_sc = OntologySinglecellLibraryConstruction()
        self.assay_type_differentiation = OntologyList(terms=["guided", "unguided"])
        self.cell_line = OntologyCellosaurus()
        self.cellontology_class = "v2021-02-01"
        self.cellontology_original = None
        self.developmental_stage = None
        self.doi = None
        self.ethnicity = None
        self.healthy = OntologyList(terms=[True, False])
        self.id = None
        self.normalization = None
        self.organ = OntologyUberon()
        self.organism = OntologyList(terms=["mouse", "human"])
        self.sample_source = OntologyList(terms=["primary_tissue", "2d_culture", "3d_culture", "cancer"])
        self.sex = OntologyList(terms=["female", "male"])
        self.year = OntologyList(terms=list(range(2000, 3000)))

    @property
    def cellontology_class(self):
        return self._cellontology_class

    @cellontology_class.setter
    def cellontology_class(self, x: str):
        self._cellontology_class = OntologyCelltypes(branch=x)
