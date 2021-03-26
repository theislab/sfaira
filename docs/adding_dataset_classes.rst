The class-based data loader python file
~~~~~~~~~~~~~~~~~~~~~~~~~~~
As an alternative to the preferred yaml-based dataloaders, users can provide a dataloader class together with the load function.
In this scenario, meta data is described in a constructor of a class in the same python file as the loading function.

1. A constructor of the following form that contains all the relevant metadata that is available before the actual dataset is loaded to memory.

.. code-block:: python

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        # Data set meta data: You do not have to include all of these and can simply skip lines corresponding
        # to attritbutes that you do not have access to. These are meta data on a sample level.
        # The meta data attributes labeled with (*) may als be supplied per cell, see below,
        # in this case, if you supply a .obs_key* attribute, you ccan leave out the sample-wise attribute.

        self.id = x  # unique identifier of data set (Organism_Organ_Year_AssaySc_NumberOfDataset_FirstAuthorLastname_doi).

        self.author = x  # author (list) who sampled / created the data set
        self.doi = x  # doi of data set accompanying manuscript

        self.download_url_data = x  # download website(s) of data files
        self.download_url_meta = x  # download website(s) of meta data files

        self.assay_sc = x  # (*, optional) protocol used to sample data (e.g. smart-seq2)
        self.assay_differentiation = x  # (*, optional) protocol used to differentiate the cell line (e.g. Lancaster, 2014)
        self.assay_type_differentiation = x  # (*, optional) type of protocol used to differentiate the cell line (guided/unguided)
        self.cell_line = x # (*, optional) cell line used (for cell culture samples)
        self.dev_stage = x  # (*, optional) developmental stage of organism
        self.ethnicity = x  # (*, optional) ethnicity of sample
        self.healthy = x  # (*, optional) whether sample represents a healthy organism
        self.normalisation = x  # (optional) normalisation applied to raw data loaded (ideally counts, "raw")
        self.organ = x  # (*, optional) organ (anatomical structure)
        self.organism = x  # (*) species / organism
        self.sample_source = x  # (*) whether the sample came from primary tissue or cell culture
        self.sex = x  # (*, optional) sex
        self.state_exact = x  # (*, optional) exact disease, treatment or perturbation state of sample
        self.year = x  # year in which sample was acquired

        # The following meta data may instead also be supplied on a cell level if an appropriate column is present in the
        # anndata instance (specifically in .obs) after loading.
        # You need to make sure this is loaded in the loading script)!
        # See above for a description what these meta data attributes mean.
        # Again, if these attributes are note available, you can simply leave this out.
        self.obs_key_assay_sc = x  # (optional, see above, do not provide if .assay_sc is provided)
        self.obs_key_assay_differentiation = x  # (optional, see above, do not provide if .age is assay_differentiation)
        self.obs_key_assay_type_differentiation = x  # (optional, see above, do not provide if .assay_type_differentiation is provided)
        self.obs_key_cell_line = x # (optional, see above, do not provide if .cell_line is provided)
        self.obs_key_dev_stage = x  # (optional, see above, do not provide if .dev_stage is provided)
        self.obs_key_ethnicity = x  # (optional, see above, do not provide if .ethnicity is provided)
        self.obs_key_healthy = x  # (optional, see above, do not provide if .healthy is provided)
        self.obs_key_organ = x  # (optional, see above, do not provide if .organ is provided)
        self.obs_key_organism = x  # (optional, see above, do not provide if .organism is provided)
        self.obs_key_sample_source = x  # (optional, see above, do not provide if .sample_source is provided)
        self.obs_key_sex = x  # (optional, see above, do not provide if .sex is provided)
        self.obs_key_state_exact = x  # (optional, see above, do not provide if .state_exact is provided)
        # Additionally, cell type annotation is ALWAYS provided per cell in .obs, this annotation is optional though.
        # name of column which contain streamlined cell ontology cell type classes:
        self.obs_key_cellontology_original = x  # (optional)
        # This cell type annotation is free text but is mapped to an ontology via a .tsv file with the same name and
        # directory as the python file of this data loader (see below).


2. A function called to load the data set into memory:
It is important to set an automated path indicating the location of the raw files here.
Our recommendation for this directory set-up is that you define a directory folder in your directory structure
in which all of these raw files will be (self.path) and then add a sub-directory named as
`self.directory_formatted_doi` (ie. the doi with all special characters replaced by "_" and place the raw files
directly into this sub directory.

.. code-block:: python

    def load(data_dir, fn=None) -> anndata.AnnData:
        fn = os.path.join(data_dir, "my.h5ad")
        adata = anndata.read(fn)  # loading instruction into adata, use other ones if the data is not h5ad
        return adata

In summary, a python file for a mouse lung data set could look like this:

.. code-block:: python

    class MyDataset(DatasetBase)
        def __init__(
                self,
                path: Union[str, None] = None,
                meta_path: Union[str, None] = None,
                cache_path: Union[str, None] = None,
                **kwargs
        ):
            super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
            self.author = "me"
            self.doi = ["my preprint", "my peer-reviewed publication"]
            self.download_url_data = "my GEO upload"
            self.normalisation = "raw"  # because I uploaded raw counts, which is good practice!
            self.organ = "lung"
            self.organism = "mouse"
            self.assay_sc = "smart-seq2"
            self.year = "2020"
            self.sample_source = "primary_tissue"

            self.obs_key_cellontology_original = "louvain_named"  # i save my cell type names in here

    def load(data_dir, fn=None) -> anndata.AnnData:
        fn = os.path.join(data_dir, "my.h5ad")
        adata = anndata.read(fn)
        return adata
