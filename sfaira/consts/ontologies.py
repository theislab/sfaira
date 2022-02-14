from typing import Dict, Union

from sfaira.versions.metadata import OntologyList, OntologyCl
from sfaira.versions.metadata import OntologyCellosaurus, OntologyHancestro, OntologyHsapdv, OntologyMondo, \
    OntologyMmusdv, OntologyEfo, OntologySex, OntologyTaxon, OntologyUberon, OntologyUberonLifecyclestage

OTHER_ORGANISM_KEY = "other"

DEFAULT_CL = "v2021-08-10"
DEFAULT_HSAPDV = "master"
DEFAULT_MONDO = "v2021-08-11"
DEFAULT_MMUSDV = "master"
DEFAULT_PATO = "v2021-08-06"
DEFAULT_NCBITAXON = "v2021-06-10"
DEFAULT_UBERON = "v2021-07-27"


class OntologyContainerSfaira:

    """
    The attributes that are relayed via properties, which have a corresponding private attribute "_*", are used to
    lazily load these ontologies upon usage and redistribute loading time from package initialisation to actual
    usage of ontology.
    """

    _assay_sc: Union[None, OntologyEfo]
    _cell_line: Union[None, OntologyCellosaurus]
    _cell_type: Union[None, OntologyCl]
    _development_stage: Union[None, Dict[str, Union[OntologyHsapdv, OntologyMmusdv]]]
    _ethnicity: Union[None, Dict[str, Union[OntologyHancestro, None]]]
    _organ: Union[None, OntologyUberon]
    _organism: Union[None, OntologyTaxon]
    _sex: Union[None, OntologySex]

    def __init__(self):
        self.key_other = OTHER_ORGANISM_KEY

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
        self.feature_type = OntologyList(terms=["rna", "peak", "protein"])
        self.id = None
        self.individual = None
        self.normalization = None
        self._organ = None
        self._organism = None
        self.primary_data = OntologyList(terms=[True, False])
        self.sample_source = OntologyList(terms=["primary_tissue", "2d_culture", "3d_culture", "tumor"])
        self._sex = None
        self.supplier = OntologyList(terms=["cellxgene", "sfaira"])
        self.tech_sample = None
        self.title = None
        self.year = OntologyList(terms=list(range(2000, 3000)))

    def reload_ontology(self, attr):
        """
        Complex alternative to attribute-wise setters.

        :param attr:
        :return:
        """
        kwargs = {"recache": True}
        if attr == "assay_sc":
            self._assay_sc = OntologyEfo(**kwargs)
        elif attr == "cell_line":
            self._cell_line = OntologyCellosaurus(**kwargs)
        elif attr == "cell_type":
            self._cell_type = OntologyCl(branch=DEFAULT_CL, **kwargs)
        elif attr == "development_stage":
            self._development_stage = {
                "Homo sapiens": OntologyHsapdv(branch=DEFAULT_HSAPDV, **kwargs),
                "Mus musculus": OntologyMmusdv(branch=DEFAULT_MMUSDV, **kwargs),
                self.key_other: OntologyUberonLifecyclestage(branch=DEFAULT_UBERON, **kwargs),
            }
        elif attr == "disease":
            self._disease = OntologyMondo(branch=DEFAULT_MONDO, **kwargs)
        elif attr == "ethnicity":
            self._ethnicity = {
                "homosapiens": OntologyHancestro(),
                self.key_other: None,
            }
        elif attr == "organ":
            self._organ = OntologyUberon(branch=DEFAULT_UBERON, **kwargs)
        elif attr == "organism":
            self._organism = OntologyTaxon(branch=DEFAULT_NCBITAXON, **kwargs)
        elif attr == "sex":
            self._sex = OntologySex(branch=DEFAULT_PATO, **kwargs)

    @property
    def assay_sc(self):
        if self._assay_sc is None:  # Lazy loading after class instantiation.
            self._assay_sc = OntologyEfo()
        return self._assay_sc

    @property
    def cell_line(self):
        if self._cell_line is None:  # Lazy loading after class instantiation.
            self._cell_line = OntologyCellosaurus()
        return self._cell_line

    @property
    def cell_type(self):
        if self._cell_type is None:  # Lazy loading after class instantiation.
            self._cell_type = OntologyCl(branch=DEFAULT_CL)
        return self._cell_type

    @cell_type.setter
    def cell_type(self, x: str):
        self._cell_type = OntologyCl(branch=x)

    @property
    def development_stage(self):
        if self._development_stage is None:  # Lazy loading after class instantiation.
            self._development_stage = {
                "Homo sapiens": OntologyHsapdv(branch=DEFAULT_HSAPDV),
                "Mus musculus": OntologyMmusdv(branch=DEFAULT_MMUSDV),
                self.key_other: OntologyUberonLifecyclestage(branch=DEFAULT_UBERON),
            }
        return self._development_stage

    @property
    def disease(self):
        if self._disease is None:  # Lazy loading after class instantiation.
            self._disease = OntologyMondo(branch=DEFAULT_MONDO)
        return self._disease

    @property
    def ethnicity(self):
        if self._ethnicity is None:  # Lazy loading after class instantiation.
            self._ethnicity = {
                "Homo sapiens": OntologyHancestro(),
                self.key_other: None,
            }
        return self._ethnicity

    @property
    def organ(self):
        if self._organ is None:  # Lazy loading after class instantiation.
            self._organ = OntologyUberon(branch=DEFAULT_UBERON)
        return self._organ

    @property
    def organism(self):
        if self._organism is None:  # Lazy loading after class instantiation.
            self._organism = OntologyTaxon(branch=DEFAULT_NCBITAXON)
        return self._organism

    @property
    def sex(self):
        if self._sex is None:  # Lazy loading after class instantiation.
            self._sex = OntologySex(branch=DEFAULT_PATO)
        return self._sex
