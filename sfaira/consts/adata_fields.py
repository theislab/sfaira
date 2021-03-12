from typing import List

"""
The classes in this file are containers of field names and element entries that are used in streamlined adata objects
in sfaira and in associated data bases.
"""


class AdataIds:
    """
    Base class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """
    age: str
    annotated: str
    assay_sc: str
    author: str
    cell_types_original: str
    cell_ontology_class: str
    cell_ontology_id: str
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
    gene_id_names: str
    healthy: str
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

    obs_keys: List[str]
    var_keys: List[str]
    uns_keys: List[str]


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
        self.cell_ontology_class = "cell_ontology_class"
        self.cell_ontology_id = "cell_ontology_id"
        self.disease = "disease"
        self.doi = "doi"
        self.dataset = "dataset"
        self.dataset_group = "dataset_group"
        self.download_url_data = "download_url_data"
        self.download_url_meta = "download_url_meta"
        self.gene_id_ensembl = "ensembl"
        self.gene_id_index = "ensembl"
        self.gene_id_names = "names"
        self.healthy = "healthy"
        self.id = "id"
        self.individual = "individual"
        self.ncells = "ncells"
        self.normalization = "normalization"
        self.organ = "organ"
        self.organism = "organism"
        self.sample_source = "sample_source"
        self.tech_sample = "tech_sample"
        self.year = "year"

        self.age = "age"
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

        self.obs_keys = [
            "age",
            "assay_sc",
            "assay_differentiation",
            "assay_type_differentiation",
            "bio_sample",
            "cell_line",
            "cell_types_original",
            "cell_ontology_class",
            "cell_ontology_id",
            "development_stage",
            "ethnicity",
            "healthy",
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
            "gene_id_names",
        ]
        self.uns_keys = [
            "annotated",
            "author",
            "doi",
            "download_url_data",
            "download_url_meta",
            "id",
            "mapped_features",
            "normalization",
            "year",
        ]


class AdataIdsCellxgene(AdataIds):
    """
    Class of constant field names of anndata.AnnData object entries", such as .uns keys and .obs columns in cellxgene
    objects.
    """
    disease_state_healthy: str
    accepted_file_names: List[str]

    def __init__(self):
        self.assay_sc = "assay"
        self.cell_types_original = "free_annotation"
        self.cell_ontology_class = "cell_type"
        self.cell_ontology_id = "cell_type_ontology_term_id"
        self.doi = "doi"
        self.disease = "disease"
        self.gene_id_names = "names"
        self.id = "id"
        self.ncells = "ncells"
        self.normalization = ""  # is always "raw"
        self.organ = "organ"
        self.organism = "organism"
        self.year = "year"

        self.age = "age"
        self.author = "contributors"
        self.development_stage = "development_stage"
        self.ethnicity = "ethnicity"
        self.sex = "sex"
        self.state_exact = "disease"

        # selected element entries used for parsing:
        self.disease_state_healthy = "normal"
        self.author_names = "names"

        # accepted file names
        self.accepted_file_names = [
            "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remixed.h5ad",
        ]

        self.obs_keys = [
            "age",
            "development_stage",
            "disease",
            "ethnicity",
            "healthy",
            "individual",
            "organ",
            "organism",
            "sex",
            "tech_sample",
        ]
        self.var_keys = [
            "gene_id_names",
        ]
        self.uns_keys = [
            "author",
            "doi",
            "id",
            "normalization",
            "year",
        ]
