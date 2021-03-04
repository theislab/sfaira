from sfaira.versions.metadata import OntologyList, OntologyCelltypes
from sfaira.versions.metadata import OntologyUberon, OntologyHsapdv, OntologyMmusdv, \
    OntologySinglecellLibraryConstruction


class OntologyContainerSfaira:

    _cellontology_class: OntologyCelltypes

    def __init__(self):
        self.ontology_age = None
        self.assay = OntologySinglecellLibraryConstruction()
        self.cellontology_class = "v2021-02-01"
        self.cellontology_original = None
        self.developmental_stage = None
        self.ethnicity = None
        self.healthy = [True, False]
        self.normalization = None
        self.organ = OntologyUberon()
        self.organism = OntologyList(terms=["mouse", "human"])
        self.sex = OntologyList(terms=["female", "male"])
        self.year = list(range(2000, 3000))

    @property
    def cellontology_class(self):
        return self._cellontology_class

    @cellontology_class.setter
    def cellontology_class(self, x: str):
        self._cellontology_class = OntologyCelltypes(branch=x)
