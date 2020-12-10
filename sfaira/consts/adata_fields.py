"""
The classes in this file are containers of field names and element entries that are used in streamlined adata objects
in sfaira and in associated data bases.
"""


class ADATA_IDS_BASE:
    """
    Base class of minimal constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """
    _annotated: str
    _author: str
    _cell_types_original: str
    _cell_ontology_class: str
    _cell_ontology_id: str
    _doi: str
    _download: str
    _dataset: str
    _dataset_group: str
    _gene_id_ensembl: str
    _gene_id_names: str
    _healthy: str
    _id: str
    _ncells: str
    _normalization: str
    _organ: str
    _protocol: str
    _species: str
    _subtissue: str
    _year: str

    @property
    def annotated(self):
        return self._annotated

    @property
    def author(self):
        return self._author

    @property
    def cell_types_original(self):
        return self._cell_types_original

    @property
    def cell_ontology_class(self):
        return self._cell_ontology_class

    @property
    def cell_ontology_id(self):
        return self._cell_ontology_id

    @property
    def dataset(self):
        return self._dataset

    @property
    def dataset_group(self):
        return self._dataset_group

    @property
    def doi(self):
        return self._doi

    @property
    def download(self):
        return self._download

    @property
    def gene_id_ensembl(self):
        return self._gene_id_ensembl

    @property
    def gene_id_names(self):
        return self._gene_id_names

    @property
    def healthy(self):
        return self._healthy

    @property
    def id(self):
        return self._id

    @property
    def ncells(self):
        return self._ncells

    @property
    def normalization(self):
        return self._normalization

    @property
    def protocol(self):
        return self._protocol

    @property
    def organ(self):
        return self._organ

    @property
    def species(self):
        return self._species

    @property
    def subtissue(self):
        return self._subtissue

    @property
    def year(self):
        return self._year


class ADATA_IDS_EXTENDED(ADATA_IDS_BASE):
    """
    Base class with extended set of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """
    _age: str
    _dev_stage: str
    _ethnicity: str
    _sex: str
    _state_exact: str

    @property
    def age(self):
        return self._age

    @property
    def dev_stage(self):
        return self._dev_stage

    @property
    def ethnicity(self):
        return self._ethnicity

    @property
    def sex(self):
        return self._sex

    @property
    def state_exact(self):
        return self._state_exact


class ADATA_IDS_SFAIRA(ADATA_IDS_EXTENDED):
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """

    def __init__(self):
        self._animal = "animal"
        self._cell_types_original = "cell_types_original"
        self._cell_ontology_class = "cell_ontology_class"
        self._cell_ontology_id = "cell_ontology_id"
        self._doi = "doi"
        self._dataset = "dataset"
        self._dataset_group = "dataset_group"
        self._gene_id_ensembl = "ensembl"
        self._gene_id_names = "names"
        self._has_celltypes = "has_celltypes"
        self._healthy = "healthy"
        self._id = "id"
        self._ncells = "ncells"
        self._normalization = "normalization"
        self._lab = "lab"
        self._organ = "organ"
        self._protocol = "protocol"
        self._subtissue = "subtissue"
        self._wget_download = "wget_download"
        self._year = "year"

        self._age = "age"
        self._dev_stage = "dev_stage"
        self._ethnicity = "ethnicity"
        self._sex = "sex"
        self._state_exact = "state_exact"


class ADATA_IDS_CELLXGENE(ADATA_IDS_EXTENDED):
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns in cellxgene
    objects.
    """
    _author_names: str
    _disease_state_healthy: str

    def __init__(self):
        self._animal = "organism"
        self._cell_types_original = "free_annotation"
        self._cell_ontology_class = "cell_type"
        self._cell_ontology_id = "cell_type_ontology_term_id"
        self._doi = ""  # TODO
        self._dataset = "dataset"
        self._dataset_group = "dataset_group"
        self._gene_id_ensembl = ""  # TODO
        self._gene_id_names = ""  # TODO
        self._has_celltypes = ""  # TODO
        self._healthy = None  # is inferred from _disease
        self._id = ""  # TODO
        self._ncells = "ncells"
        self._normalization = None  # is always "counts"
        self._lab = ""  # TODO
        self._organ = ""  # TODO
        self._protocol = "assay"
        self._subtissue = ""  # TODO
        self._wget_download = ""  # TODO
        self._year = ""  # TODO

        self._age = "age"
        self._author = "contributors"
        self._dev_stage = "development_stage"
        self._ethnicity = "ethnicity"
        self._sex = "sex"
        self._state_exact = "disease"

        # selected element entries used for parsing:
        self._disease_state_healthy = "normal"
        self._author_names = "names"

    @property
    def author_names(self):
        return self._author_names

    @property
    def disease_state_healthy(self):
        return self._disease_state_healthy
