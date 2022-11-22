from typing import Dict, Union

from sfaira.consts.schema import DEFAULT_SCHEMA, ONTOLOGY_VERSIONS
from sfaira.versions.metadata import OntologyList, OntologyCl
from sfaira.versions.metadata import OntologyCellosaurus, OntologyHancestro, OntologyHsapdv, OntologyMondo, \
    OntologyMmusdv, OntologyEfo, OntologyPato, OntologyTaxon, OntologyUberon, OntologyUberonLifecyclestage

OTHER_ORGANISM_KEY = "other"


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
    _sex: Union[None, OntologyPato]

    versions: dict

    def __init__(self):
        self.key_other = OTHER_ORGANISM_KEY
        self.set_schema_version()

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
        self._organ = None
        self._organism = None
        self.primary_data = OntologyList(terms=[True, False])
        self.sample_source = OntologyList(terms=["primary_tissue", "2d_culture", "3d_culture", "tumor"])
        self.suspension_type = OntologyList(terms=["cell", "na", "nucleus"])
        self._sex = None
        self.supplier = OntologyList(terms=["cellxgene", "sfaira"])
        self.tech_sample = None
        self.title = None
        self.year = OntologyList(terms=list(range(2000, 3000)))

    @property
    def versioned_ontologies(self):
        return [
            "assay_sc",
            "cell_line",
            "cell_type",
            "development_stage",
            "disease",
            "ethnicity",
            "organ",
            "organism",
            "sex"
        ]

    def set_schema_version(self, version: str = DEFAULT_SCHEMA):
        self.versions = ONTOLOGY_VERSIONS[version]
        for k in self.versioned_ontologies:
            self.reload_ontology(attr=k)

    def reload_ontology(self, attr, recache: bool = False):
        """
        Complex alternative to attribute-wise setters.

        :param attr:
        :param recache: Whether to re-download and cache the ontology.
        :return:
        """
        kwargs = {"recache": recache}
        if attr == "assay_sc":
            self._assay_sc = OntologyEfo(branch=self.versions["VERSION_EFO"], **kwargs)
        elif attr == "cell_line":
            self._cell_line = OntologyCellosaurus(**kwargs)
        elif attr == "cell_type":
            self._cell_type = OntologyCl(branch=self.versions["VERSION_CL"], **kwargs)
        elif attr == "development_stage":
            self._development_stage = {
                "Homo sapiens": OntologyHsapdv(branch=self.versions["VERSION_HSAPDV"], **kwargs),
                "Mus musculus": OntologyMmusdv(branch=self.versions["VERSION_MMUSDV"], **kwargs),
                self.key_other: OntologyUberonLifecyclestage(branch=self.versions["VERSION_UBERON"], **kwargs),
            }
        elif attr == "disease":
            self._disease = OntologyMondo(branch=self.versions["VERSION_MONDO"], **kwargs)
        elif attr == "ethnicity":
            self._ethnicity = {
                "Homo sapiens": OntologyHancestro(branch=self.versions["VERSION_HANCESTRO"]),
                self.key_other: None,
            }
        elif attr == "organ":
            self._organ = OntologyUberon(branch=self.versions["VERSION_UBERON"], **kwargs)
        elif attr == "organism":
            self._organism = OntologyTaxon(branch=self.versions["VERSION_NCBITAXON"], **kwargs)
        elif attr == "sex":
            self._sex = OntologyPato(branch=self.versions["VERSION_PATO"], **kwargs)

    @property
    def assay_sc(self):
        if self._assay_sc is None:  # Lazy loading after class instantiation.
            self._assay_sc = OntologyEfo(branch=self.versions["VERSION_EFO"])
        return self._assay_sc

    @assay_sc.setter
    def assay_sc(self, x: str):
        self._assay_sc = OntologyEfo(branch=x)

    @property
    def cell_line(self):
        if self._cell_line is None:  # Lazy loading after class instantiation.
            self._cell_line = OntologyCellosaurus()
        return self._cell_line

    @property
    def cell_type(self):
        if self._cell_type is None:  # Lazy loading after class instantiation.
            self._cell_type = OntologyCl(branch=self.versions["VERSION_CL"])
        return self._cell_type

    @cell_type.setter
    def cell_type(self, x: str):
        self._cell_type = OntologyCl(branch=x)

    @property
    def development_stage(self):
        if self._development_stage is None:  # Lazy loading after class instantiation.
            self._development_stage = {
                "Homo sapiens": OntologyHsapdv(branch=self.versions["VERSION_HSAPDV"]),
                "Mus musculus": OntologyMmusdv(branch=self.versions["VERSION_MMUSDV"]),
                self.key_other: OntologyUberonLifecyclestage(branch=self.versions["VERSION_UBERON"]),
            }
        return self._development_stage

    @property
    def disease(self):
        if self._disease is None:  # Lazy loading after class instantiation.
            self._disease = OntologyMondo(branch=self.versions["VERSION_MONDO"])
        return self._disease

    @disease.setter
    def disease(self, x: str):
        self._disease = OntologyMondo(branch=x)

    @property
    def ethnicity(self):
        if self._ethnicity is None:  # Lazy loading after class instantiation.
            self._ethnicity = {
                "Homo sapiens": OntologyHancestro(branch=self.versions["VERSION_HANCESTRO"]),
                self.key_other: None,
            }
        return self._ethnicity

    @property
    def organ(self):
        if self._organ is None:  # Lazy loading after class instantiation.
            self._organ = OntologyUberon(branch=self.versions["VERSION_UBERON"])
        return self._organ

    @organ.setter
    def organ(self, x: str):
        self._organ = OntologyUberon(branch=x)

    @property
    def organism(self):
        if self._organism is None:  # Lazy loading after class instantiation.
            self._organism = OntologyTaxon(branch=self.versions["VERSION_NCBITAXON"])
        return self._organism

    @organism.setter
    def organism(self, x: str):
        self._organism = OntologyTaxon(branch=x)

    @property
    def sex(self):
        if self._sex is None:  # Lazy loading after class instantiation.
            self._sex = OntologyPato(branch=self.versions["VERSION_PATO"])
        return self._sex

    @sex.setter
    def sex(self, x: str):
        self._sex = OntologyPato(branch=x)


OC = OntologyContainerSfaira()
