from sfaira.versions.metadata import OntologyList, OntologyCl
from sfaira.versions.metadata import OntologyCellosaurus, OntologyHsapdv, OntologyMondo, \
    OntologyMmusdv, OntologySinglecellLibraryConstruction, OntologyUberon


class OntologyContainerSfaira:

    _cellontology_class: OntologyCl

    def __init__(self):
        self.annotated = OntologyList(terms=[True, False])
        self.author = None
        self.assay_differentiation = None
        self.assay_sc = OntologySinglecellLibraryConstruction()
        self.assay_type_differentiation = OntologyList(terms=["guided", "unguided"])
        self.bio_sample = None
        self.cell_line = OntologyCellosaurus()
        self.cellontology_class = "v2021-02-01"
        self.cell_types_original = None
        self.collection_id = None
        self.default_embedding = None
        self.development_stage = {
            "human": OntologyHsapdv(),
            "mouse": OntologyMmusdv(),
        }
        self.disease = OntologyMondo()
        self.doi = None
        self.ethnicity = {
            "human": None,  # TODO OntologyHancestro
            "mouse": None,
        }
        self.id = None
        self.individual = None
        self.normalization = None
        self.organ = OntologyUberon()
        self.organism = OntologyList(terms=["mouse", "human"])  # TODO introduce NCBItaxon here
        self.primary_data = OntologyList(terms=[True, False])
        self.sample_source = OntologyList(terms=["primary_tissue", "2d_culture", "3d_culture", "tumor"])
        self.sex = OntologyList(terms=["female", "male", "mixed", "unknown", "other"])
        self.supplier = OntologyList(terms=["cellxgene", "sfaira"])
        self.tech_sample = None
        self.title = None
        self.year = OntologyList(terms=list(range(2000, 3000)))

    @property
    def cellontology_class(self):
        return self._cellontology_class

    @cellontology_class.setter
    def cellontology_class(self, x: str):
        self._cellontology_class = OntologyCl(branch=x)
