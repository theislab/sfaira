from typing import Dict, Union

from sfaira.versions.metadata import OntologyList, OntologyCl
from sfaira.versions.metadata import OntologyCellosaurus, OntologyHancestro, OntologyHsapdv, OntologyMondo, \
    OntologyMmusdv, OntologySinglecellLibraryConstruction, OntologyUberon

DEFAULT_CL = "v2021-02-01"
DEFAULT_UBERON = "2019-11-22"


class OntologyContainerSfaira:

    """
    The attributes that are relayed via properties, which have a corresponding private attribute "_*", are used to
    lazily load these ontologies upon usage and redistribute loading time from package initialisation to actual
    usage of ontology.
    """

    _assay_sc: Union[None, OntologySinglecellLibraryConstruction]
    _cell_line: Union[None, OntologyCellosaurus]
    _cell_type: Union[None, OntologyCl]
    _development_stage: Union[None, Dict[str, Union[OntologyHsapdv, OntologyMmusdv]]]
    _ethnicity: Union[None, Dict[str, Union[OntologyHancestro, None]]]
    _organ: Union[None, OntologyUberon]

    def __init__(self):
        self.annotated = OntologyList(terms=[True, False])
        self.author = None
        self.assay_differentiation = None
        self._assay_sc = None
        self.assay_type_differentiation = OntologyList(terms=["guided", "unguided"])
        self.bio_sample = None
        self._cell_line = None
        self._cell_type = None
        self.collection_id = None
        self.default_embedding = None
        self._development_stage = None
        self._disease = None
        self.doi = None
        self.doi_main = None
        self.doi_journal = None
        self.doi_preprint = None
        self._ethnicity = None
        self.id = None
        self.individual = None
        self.normalization = None
        self._organ = None
        self.organism = OntologyList(terms=["mouse", "human"])  # TODO introduce NCBItaxon here
        self.primary_data = OntologyList(terms=[True, False])
        self.sample_source = OntologyList(terms=["primary_tissue", "2d_culture", "3d_culture", "tumor"])
        self.sex = OntologyList(terms=["female", "male", "mixed", "unknown", "other"])
        self.supplier = OntologyList(terms=["cellxgene", "sfaira"])
        self.tech_sample = None
        self.title = None
        self.year = OntologyList(terms=list(range(2000, 3000)))

    def reload_ontology(self, attr):
        kwargs = {"recache": True}
        if attr == "assay_sc":
            self._assay_sc = OntologySinglecellLibraryConstruction(**kwargs)
        elif attr == "cell_line":
            self._cell_line = OntologyCellosaurus(**kwargs)
        elif attr == "cellontology_class":
            self._cell_type = OntologyCl(branch=DEFAULT_CL, **kwargs)
        elif attr == "development_stage":
            self._development_stage = {
                "human": OntologyHsapdv(),
                "mouse": OntologyMmusdv(),
            }
        elif attr == "disease":
            self._disease = OntologyMondo(**kwargs)
        elif attr == "ethnicity":
            self._ethnicity = {
                "human": OntologyHancestro(),
                "mouse": None,
            }
        elif attr == "organ":
            self._organ = OntologyUberon(**kwargs)
        return self._assay_sc

    @property
    def assay_sc(self):
        if self._assay_sc is None:
            self._assay_sc = OntologySinglecellLibraryConstruction()
        return self._assay_sc

    @property
    def cell_line(self):
        if self._cell_line is None:
            self._cell_line = OntologyCellosaurus()
        return self._cell_line

    @property
    def cell_type(self):
        if self._cell_type is None:
            self._cell_type = OntologyCl(branch=DEFAULT_CL)
        return self._cell_type

    @cell_type.setter
    def cell_type(self, x: str):
        self._cell_type = OntologyCl(branch=x)

    @property
    def development_stage(self):
        if self._development_stage is None:
            self._development_stage = {
                "human": OntologyHsapdv(),
                "mouse": OntologyMmusdv(),
            }
        return self._development_stage

    @property
    def disease(self):
        if self._disease is None:
            self._disease = OntologyMondo()
        return self._disease

    @property
    def ethnicity(self):
        if self._ethnicity is None:
            self._ethnicity = {
                "human": OntologyHancestro(),
                "mouse": None,
            }
        return self._ethnicity

    @property
    def organ(self):
        if self._organ is None:
            self._organ = OntologyUberon(branch=DEFAULT_UBERON)
        return self._organ
