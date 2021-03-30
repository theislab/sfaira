from sfaira.versions.metadata import OntologyList, OntologyCelltypes
from sfaira.versions.metadata import OntologyUberon, OntologyHsapdv, OntologyMondo, OntologyMmusdv, \
    OntologySinglecellLibraryConstruction, OntologyCellosaurus


class OntologyContainerSfaira:

    _cellontology_class: OntologyCelltypes

    def __init__(self):
        self.annotated = OntologyList(terms=[True, False])
        self.author = None
        self.assay_differentiation = None
        self.assay_sc = OntologySinglecellLibraryConstruction()
        self.assay_type_differentiation = OntologyList(terms=["guided", "unguided"])
        self.cell_line = OntologyCellosaurus()
        self.cellontology_class = "v2021-02-01"
        self.cellontology_original = None
        self.development_stage = None  # OntologyHsapdv()  # TODO allow for other organisms here too.
        self.disease = OntologyMondo()
        self.doi = None
        self.ethnicity = None  # OntologyHancestro()
        self.healthy = OntologyList(terms=[True, False])
        self.id = None
        self.normalization = None
        self.organ = OntologyUberon()
        self.organism = OntologyList(terms=["mouse", "human"])  # TODO introduce NCBItaxon here
        self.sample_source = OntologyList(terms=["primary_tissue", "2d_culture", "3d_culture", "tumor"])
        self.sex = OntologyList(terms=["female", "male", "mixed", "unknown", "other"])
        self.year = OntologyList(terms=list(range(2000, 3000)))

    @property
    def cellontology_class(self):
        return self._cellontology_class

    @cellontology_class.setter
    def cellontology_class(self, x: str):
        self._cellontology_class = OntologyCelltypes(branch=x)
