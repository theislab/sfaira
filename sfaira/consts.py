
class ADATA_IDS:
    """
    Class of constant field names of anndata.AnnData object entries, such as .uns keys and .obs columns.
    """

    def __init__(self):
        self._age = "age"
        self._animal = "animal"
        self._cell_types_original = "cell_types_original"
        self._cell_ontology_class = "cell_ontology_class"
        self._cell_ontology_id = "cell_ontology_id"
        self._dev_stage = "dev_stage"
        self._doi = "doi"
        self._dataset = "dataset"
        self._dataset_group = "dataset_group"
        self._gene_id_ensembl = "ensembl"
        self._gene_id_names = "names"
        self._has_celltypes = "has_celltypes"
        self._healthy = "healthy"
        self._id = "id"
        self._normalization = "normalization"
        self._lab = "lab"
        self._organ = "organ"
        self._protocol = "protocol"
        self._state_exact = "state_exact"
        self._subtissue = "subtissue"
        self._wget_download = "wget_download"
        self._year = "year"

    @property
    def age(self):
        return self._age

    @property
    def animal(self):
        return self._animal

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
    def dev_stage(self):
        return self._dev_stage

    @property
    def doi(self):
        return self._doi

    @property
    def gene_id_ensembl(self):
        return self._gene_id_ensembl

    @property
    def gene_id_names(self):
        return self._gene_id_names

    @property
    def has_celltypes(self):
        return self._has_celltypes

    @property
    def healthy(self):
        return self._healthy

    @property
    def id(self):
        return self._id

    @property
    def lab(self):
        return self._lab

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
    def subtissue(self):
        return self._subtissue

    @property
    def state_exact(self):
        return self._state_exact

    @property
    def wget_download(self):
        return self._wget_download

    @property
    def year(self):
        return self._year

