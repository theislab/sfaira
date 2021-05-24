"""
The classes in this file are containers of field names and element entries that are used in streamlined adata objects
in sfaira and in associated data bases.
"""
from typing import List, Union


class AdataIds:
    """
    Base class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """
    annotated: str
    assay_sc: str
    author: str
    cell_types_original: str
    cellontology_class: str
    cellontology_id: str
    development_stage: str
    disease: str
    doi: str
    download_url_data: str
    download_url_meta: str
    dataset: str
    dataset_group: str
    ethnicity: str
    gene_id_ensembl: str
    gene_id_index: str
    gene_id_symbols: str
    id: str
    individual: str
    ncells: str
    normalization: str
    organ: str
    organism: str
    sample_source: str
    sex: str
    state_exact: str
    tech_sample: str
    year: str

    load_raw: str
    mapped_features: str
    remove_gene_version: str

    obs_keys: List[str]
    var_keys: List[str]
    uns_keys: List[str]

    classmap_source_key: str
    classmap_target_key: str
    classmap_target_id_key: str

    unknown_celltype_identifier: Union[str, None]
    not_a_cell_celltype_identifier: Union[str, None]
    unknown_metadata_identifier: Union[str, None]

    @property
    def controlled_meta_keys(self):
        return [getattr(self, k) for k in self.obs_keys + self.uns_keys]

    @property
    def controlled_meta_fields(self):
        return [k for k in self.obs_keys + self.uns_keys]


class AdataIdsSfaira(AdataIds):
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns in sfaira
    dataloader objects.
    """

    assay_differentiation: str
    assay_type_differentiation: str
    bio_sample: str
    cell_line: str

    def __init__(self):
        self.annotated = "annotated"
        self.assay_sc = "assay_sc"
        self.assay_differentiation = "assay_differentiation"
        self.assay_type_differentiation = "assay_type_differentiation"
        self.author = "author"
        self.bio_sample = "bio_sample"
        self.cell_line = "cell_line"
        self.cell_types_original = "cell_types_original"
        self.cellontology_class = "cell_ontology_class"
        self.cellontology_id = "cell_ontology_id"
        self.default_embedding = "default_embedding"
        self.disease = "disease"
        self.doi = "doi"
        self.dataset = "dataset"
        self.dataset_group = "dataset_group"
        self.download_url_data = "download_url_data"
        self.download_url_meta = "download_url_meta"
        self.gene_id_ensembl = "ensembl"
        self.gene_id_index = self.gene_id_ensembl
        self.gene_id_symbols = "names"
        self.id = "id"
        self.individual = "individual"
        self.ncells = "ncells"
        self.normalization = "normalization"
        self.organ = "organ"
        self.organism = "organism"
        self.primary_data = "primary_data"
        self.sample_source = "sample_source"
        self.tech_sample = "tech_sample"
        self.title = "title"
        self.year = "year"

        self.development_stage = "development_stage"
        self.ethnicity = "ethnicity"
        self.sex = "sex"
        self.state_exact = "state_exact"

        self.load_raw = "load_raw"
        self.mapped_features = "mapped_features"
        self.remove_gene_version = "remove_gene_version"

        self.classmap_source_key = "source"
        self.classmap_target_key = "target"
        self.classmap_target_id_key = "target_id"

        self.unknown_celltype_identifier = "UNKNOWN"
        self.not_a_cell_celltype_identifier = "NOT_A_CELL"
        self.unknown_metadata_identifier = "unknown"

        self.obs_keys = [
            "assay_sc",
            "assay_differentiation",
            "assay_type_differentiation",
            "bio_sample",
            "cell_line",
            "cell_types_original",
            "cellontology_class",
            "cellontology_id",
            "development_stage",
            "disease",
            "ethnicity",
            "individual",
            "organ",
            "organism",
            "sex",
            "state_exact",
            "sample_source",
            "tech_sample",
        ]
        self.var_keys = [
            "gene_id_ensembl",
            "gene_id_symbols",
        ]
        self.uns_keys = [
            "annotated",
            "author",
            "default_embedding",
            "doi",
            "download_url_data",
            "download_url_meta",
            "id",
            "mapped_features",
            "ncells",
            "normalization",
            "primary_data",
            "title",
            "year",
            "load_raw",
            "mapped_features",
            "remove_gene_version",
        ]


class AdataIdsCellxgene(AdataIds):
    """
    Class of constant field names of anndata.AnnData object entries", such as .uns keys and .obs columns in cellxgene
    objects.
    """
    accepted_file_names: List[str]

    def __init__(self):
        self.assay_sc = "assay"
        self.cell_types_original = "free_annotation"  # TODO "free_annotation" not always given
        # TODO: -> This will break streamlining though if self.cell_types_original is the same value as self.cellontology_class!!
        self.cellontology_class = "cell_type"
        self.cellontology_id = "cell_type_ontology_term_id"
        self.default_embedding = "default_embedding"
        self.doi = "preprint_doi"
        self.disease = "disease"
        self.gene_id_symbols = "gene_symbol"
        self.gene_id_index = self.gene_id_symbols
        self.id = "id"
        self.ncells = "ncells"
        self.organ = "tissue"
        self.organism = "organism"
        self.title = "title"
        self.year = "year"

        self.author = "contributors"
        self.development_stage = "development_stage"
        self.ethnicity = "ethnicity"
        self.sex = "sex"
        self.state_exact = "disease"
        self.tech_sample = "batch"

        # selected element entries used for parsing:
        self.author_names = "names"

        self.unknown_celltype_identifier = None
        self.not_a_cell_celltype_identifier = self.unknown_celltype_identifier
        self.unknown_metadata_identifier = "unknown"
        self.invalid_metadata_identifier = "na"
        self.unknown_metadata_ontology_id_identifier = ""

        # accepted file names
        self.accepted_file_names = [
            "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remixed.h5ad",
        ]

        self.obs_keys = [
            "assay_sc",
            "cell_types_original",
            "cellontology_class",
            "cellontology_id",
            "development_stage",
            "disease",
            "ethnicity",
            "organ",
            "organism",
            "sex",
            "tech_sample",
        ]
        self.var_keys = [
            "gene_id_symbols",
        ]
        self.uns_keys = [
            "default_embedding",
            "id",
            "title",
        ]
        # These attributes related to obs and uns keys above are also in the data set attributes that can be
        # inquired before download via the REST API:
        self.dataset_keys = [
            "assay_sc",
            "development_stage",
            "disease",
            "ethnicity",
            "organ",
            "organism",
            "sex",
        ]
