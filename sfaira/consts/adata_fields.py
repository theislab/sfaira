from typing import List

"""
The classes in this file are containers of field names and element entries that are used in streamlined adata objects
in sfaira and in associated data bases.
"""


class AdataIds:
    """
    Base class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """
    _age: str
    _annotated: str
    _assay_sc: str
    _assay_differentiation: str
    _assay_type_differentiation: str
    _author: str
    _bio_sample: str
    _cell_line: str
    _cell_types_original: str
    _cell_ontology_class: str
    _cell_ontology_id: str
    _development_stage: str
    _doi: str
    _download_url_data: str
    _download_url_meta: str
    _dataset: str
    _dataset_group: str
    _ethnicity: str
    _gene_id_ensembl: str
    _gene_id_index: str
    _gene_id_names: str
    _healthy: str
    _id: str
    _individual: str
    _ncells: str
    _normalization: str
    _organ: str
    _organism: str
    _sample_source: str
    _sex: str
    _state_exact: str
    _tech_sample: str
    _year: str

    @property
    def age(self) -> str:
        return self._age

    @property
    def annotated(self) -> str:
        return self._annotated

    @property
    def assay_sc(self) -> str:
        return self._assay_sc

    @property
    def assay_differentiation(self) -> str:
        return self._assay_differentiation

    @property
    def assay_type_differentiation(self) -> str:
        return self._assay_type_differentiation

    @property
    def author(self) -> str:
        return self._author

    @property
    def bio_sample(self) -> str:
        return self._bio_sample

    @property
    def cell_line(self) -> str:
        return self._cell_line

    @property
    def cell_types_original(self) -> str:
        return self._cell_types_original

    @property
    def cell_ontology_class(self) -> str:
        return self._cell_ontology_class

    @property
    def cell_ontology_id(self) -> str:
        return self._cell_ontology_id

    @property
    def dataset(self) -> str:
        return self._dataset

    @property
    def dataset_group(self) -> str:
        return self._dataset_group

    @property
    def development_stage(self) -> str:
        return self._development_stage

    @property
    def doi(self) -> str:
        return self._doi

    @property
    def download_url_data(self) -> str:
        return self._download_url_data

    @property
    def download_url_meta(self) -> str:
        return self._download_url_meta

    @property
    def ethnicity(self) -> str:
        return self._ethnicity

    @property
    def gene_id_ensembl(self) -> str:
        return self._gene_id_ensembl

    @property
    def gene_id_index(self) -> str:
        return self._gene_id_index

    @gene_id_index.setter
    def gene_id_index(self, x: str):
        self._gene_id_index = x

    @property
    def gene_id_names(self) -> str:
        return self._gene_id_names

    @property
    def healthy(self) -> str:
        return self._healthy

    @property
    def id(self) -> str:
        return self._id

    @property
    def individual(self) -> str:
        return self._individual

    @property
    def ncells(self) -> str:
        return self._ncells

    @property
    def normalization(self) -> str:
        return self._normalization

    @property
    def organ(self) -> str:
        return self._organ

    @property
    def organism(self) -> str:
        return self._organism

    @property
    def sample_source(self) -> str:
        return self._sample_source

    @property
    def sex(self) -> str:
        return self._sex

    @property
    def state_exact(self) -> str:
        return self._state_exact

    @property
    def tech_sample(self) -> str:
        return self._tech_sample

    @property
    def year(self) -> str:
        return self._year


class AdataIdsSfaira(AdataIds):
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns in sfaira
    dataloader objects.
    """

    def __init__(self):
        self._annotated = "annotated"
        self._assay_sc = "assay_sc"
        self._assay_differentiation = "assay_differentiation"
        self._assay_type_differentiation = "assay_type_differentiation"
        self._author = "author"
        self._bio_sample = "bio_sample"
        self._cell_line = "cell_line"
        self._cell_types_original = "cell_types_original"
        self._cell_ontology_class = "cell_ontology_class"
        self._cell_ontology_id = "cell_ontology_id"
        self._doi = "doi"
        self._dataset = "dataset"
        self._dataset_group = "dataset_group"
        self._download_url_data = "download_url_data"
        self._download_url_meta = "download_url_meta"
        self._gene_id_ensembl = "ensembl"
        self._gene_id_index = "ensembl"
        self._gene_id_names = "names"
        self._healthy = "healthy"
        self._id = "id"
        self._individual = "individual"
        self._ncells = "ncells"
        self._normalization = "normalization"
        self._organ = "organ"
        self._organism = "organism"
        self._sample_source = "sample_source"
        self._tech_sample = "bio_sample"
        self._year = "year"

        self._age = "age"
        self._development_stage = "development_stage"
        self._ethnicity = "ethnicity"
        self._sex = "sex"
        self._state_exact = "state_exact"

        self._load_raw = "load_raw"
        self._mapped_features = "mapped_features"
        self._remove_gene_version = "remove_gene_version"

        self.classmap_source_key = "source"
        self.classmap_target_key = "target"
        self.classmap_target_id_key = "target_id"

        self.unknown_celltype_identifier = "UNKNOWN"
        self.not_a_cell_celltype_identifier = "NOT_A_CELL"

    @property
    def load_raw(self) -> str:
        return self._load_raw

    @property
    def mapped_features(self) -> str:
        return self._mapped_features

    @property
    def remove_gene_version(self) -> str:
        return self._remove_gene_version


class AdataIdsCellxgene(AdataIds):
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns in cellxgene
    objects.
    """
    _author_names: str
    _disease_state_healthy: str
    accepted_file_names: List[str]

    def __init__(self):
        self._assay_sc = "assay"
        self._cell_types_original = "free_annotation"
        self._cell_ontology_class = "cell_type"
        self._cell_ontology_id = "cell_type_ontology_term_id"
        self._doi = ""  # TODO
        self._dataset = "dataset"
        self._dataset_group = "dataset_group"
        self._download_url_data = ""  # TODO
        self._download_url_meta = ""  # never necessary as we interface via anndata objects
        self._gene_id_ensembl = ""  # TODO
        self._gene_id_index = "ensembl"
        self._gene_id_names = ""  # TODO
        self._has_celltypes = ""  # TODO
        self._healthy = None  # is inferred from _disease
        self._id = ""  # TODO
        self._ncells = "ncells"
        self._normalization = ""  # is always "raw"
        self._organ = ""  # TODO
        self._organism = "organism"
        self._year = ""  # TODO

        self._age = "age"
        self._author = "contributors"
        self._development_stage = "development_stage"
        self._ethnicity = "ethnicity"
        self._sex = "sex"
        self._state_exact = "disease"

        # selected element entries used for parsing:
        self._disease_state_healthy = "normal"
        self._author_names = "names"

        # accepted file names
        self.accepted_file_names = [
            "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remixed.h5ad",
        ]

    @property
    def author_names(self) -> str:
        return self._author_names

    @property
    def disease_state_healthy(self) -> str:
        return self._disease_state_healthy
