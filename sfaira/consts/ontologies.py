from sfaira.versions.metadata import OntologyList, OntologyCelltypes
from sfaira.versions.metadata import OntologyUberon, OntologyHsapdv, OntologyMmusdv, \
    OntologySinglecellLibraryConstruction


class OntologyContainerSfaira:

    def __init__(self):
        self.ontology_age = None
        self._ontology_cell_types = None
        self.ontology_cell_types = "v2021-02-01"
        self.ontology_dev_stage = None
        self.ontology_ethnicity = None
        self.ontology_healthy = [True, False]
        self.ontology_normalization = None
        self.ontology_organ = OntologyUberon()
        self.ontology_organism = OntologyList(terms=["mouse", "human"])
        self.ontology_protocol = OntologySinglecellLibraryConstruction()
        self.ontology_sex = OntologyList(terms=["female", "male"])
        self.ontology_subtissue = None
        self.ontology_year = list(range(2000, 3000))

    @property
    def ontology_cell_types(self):
        return self._ontology_cell_types

    @ontology_cell_types.setter
    def ontology_cell_types(self, x: str):
        self._ontology_cell_types = OntologyCelltypes(branch=x)
