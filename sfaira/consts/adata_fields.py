import numpy as np
from typing import List

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
    _download_meta: str
    _dataset: str
    _dataset_group: str
    _gene_id_ensembl: str
    _gene_id_index: str
    _gene_id_names: str
    _healthy: str
    _id: str
    _ncells: str
    _normalization: str
    _organ: str
    _organism: str
    _protocol: str
    _subtissue: str
    _year: str

    @property
    def annotated(self) -> str:
        return self._annotated

    @property
    def author(self) -> str:
        return self._author

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
    def doi(self) -> str:
        return self._doi

    @property
    def download(self) -> str:
        return self._download

    @property
    def download_meta(self) -> str:
        return self._download_meta

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
    def ncells(self) -> str:
        return self._ncells

    @property
    def normalization(self) -> str:
        return self._normalization

    @property
    def organ(self) -> str:
        return self._organ

    @property
    def organism(self) -> str:  # TODO refactor into organism
        return self._organism

    @property
    def protocol(self) -> str:
        return self._protocol

    @property
    def subtissue(self) -> str:
        return self._subtissue

    @property
    def year(self) -> str:
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
    def age(self) -> str:
        return self._age

    @property
    def dev_stage(self) -> str:
        return self._dev_stage

    @property
    def ethnicity(self) -> str:
        return self._ethnicity

    @property
    def sex(self) -> str:
        return self._sex

    @property
    def state_exact(self) -> str:
        return self._state_exact


class ADATA_IDS_SFAIRA(ADATA_IDS_EXTENDED):
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """

    def __init__(self):
        self._annotated = "annotated"
        self._author = "author"
        self._cell_types_original = "cell_types_original"
        self._cell_ontology_class = "cell_ontology_class"
        self._cell_ontology_id = "cell_ontology_id"
        self._doi = "doi"
        self._dataset = "dataset"
        self._dataset_group = "dataset_group"
        self._download = "download"
        self._download_meta = "download_meta"
        self._gene_id_ensembl = "ensembl"
        self._gene_id_index = "ensembl"
        self._gene_id_names = "names"
        self._healthy = "healthy"
        self._id = "id"
        self._ncells = "ncells"
        self._normalization = "normalization"
        self._organ = "organ"
        self._organism = "organism"
        self._protocol = "protocol"
        self._subtissue = "subtissue"
        self._year = "year"

        self._age = "age"
        self._dev_stage = "dev_stage"
        self._ethnicity = "ethnicity"
        self._sex = "sex"
        self._state_exact = "state_exact"

        self._load_raw = "load_raw"
        self._mapped_features = "mapped_features"
        self._remove_gene_version = "remove_gene_version"

        # Allowed field values:
        self.age_allowed_entries = None
        self.dev_stage_allowed_entries = None
        self.ethnicity_allowed_entries = None
        self.normalization_allowed_entries = None
        self.organ_allowed_entries = None
        self.organism_allowed_entries = ["mouse", "human"]
        self.protocol_allowed_entries = None
        self.sex_allowed_entries = None
        self.subtissue_allowed_entries = None
        self.year_allowed_entries = list(range(2000, 3000))
        # Free fields that are not constrained:
        # _author, _download, _download_meta, _doi, _id, _state_exact

        self.unknown_celltype_name = "unknown"
        self.unknown_celltype_identifiers = ["nan", "none", "unknown", np.nan, None]

    @property
    def load_raw(self) -> str:
        return self._load_raw

    @property
    def mapped_features(self) -> str:
        return self._mapped_features

    @property
    def remove_gene_version(self) -> str:
        return self._remove_gene_version


class ADATA_IDS_CELLXGENE(ADATA_IDS_EXTENDED):
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns in cellxgene
    objects.
    """
    _author_names: str
    _disease_state_healthy: str
    accepted_file_names: List[str]

    def __init__(self):
        self._cell_types_original = "free_annotation"
        self._cell_ontology_class = "cell_type"
        self._cell_ontology_id = "cell_type_ontology_term_id"
        self._doi = ""  # TODO
        self._dataset = "dataset"
        self._dataset_group = "dataset_group"
        self._download = ""  # TODO
        self._download_meta = ""  # never necessary as we interface via anndata objects
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
        self._protocol = "assay"
        self._subtissue = ""  # TODO
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
