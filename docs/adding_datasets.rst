Adding data sets
===================

Adding datasets to sfaira is a great way to increase the visibility of your dataset and to make it available to a large audience.
This process requires a couple of steps as outlined in the following sections.

    1. Write a dataloader as outlined below.
    2. Identify the raw files as indicated in the dataloader classes and copy them into your directory structure as required by your data loader.
    3. You can contribute the data loader to public sfaira, we do not manage data upload though.
       During publication, you would upload this data set to a server like GEO and the data loader contributed to sfaira would use this download link.

The following sections will first describe the underlying design principles of sfaira dataloaders and
then explain how to interactively create, validate and test dataloaders.

Use data loaders on existing data repository
--------------------------------------------

You only want to use data sets with existing data loaders and have adapted your directory structure as above?
In that case, you can immediately start using the data loader functions, you just need to supply the root directory
of the directory structure as `path to the constructor of the class that you are using.
Depending on the functionalities you want to use, you would often want to create a directory with cached meta data
first. This can be easily done via the script sfaira.data.utils.create_meta.py. This meta information is necessary to
anticipate file sizes for backing merged adata objects, for example, and is used for lazy loading.

Writing dataloaders
---------------------

The study-centric data loader module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the sfaira code, data loaders are organised into directories, which correspond to publications.
All data loaders corresponding to data sets of one study are grouped into this directory.
Next, each data set is represented by one data loader python file in this directory.
See below for more complex set ups with repetitive data loader code.

Check that the data loader was not already implemented
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will open issues for all planned data loaders, so you can search both the code_ base and our GitHub issues_ for
matching data loaders before you start writing one.
The core data loader identified is the directory compatible doi,
which is the doi with all special characters replaced by "_" and a "d" prefix is used:
"10.1016/j.cell.2019.06.029" becomes "d10_1016_j_cell_2019_06_029".
Searching for this string should yield a match if it is already implemented, take care to look for both
preprint and publication DOIs if both are available. We will also mention publication names in issues, you will however not find these in the code.

.. _code: https://github.com/theislab/sfaira/tree/dev
.. _issues: https://github.com/theislab/sfaira/issues


The data loader python file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each data set (organsism, organ, assay_sc, optionally also batches) has its own data loader class. Each such class is
in a separate file and inherits from a base class that contains most functionalities. Accordingly, the data loader class
looks very similar in parts to a cell in a juypter notebook that performs data loading. The core features that must be included are:

1. A constructor of the following form that can be used to interact with the data set
before it is loaded into memory:

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

        self.age = x  # (*, optional) age of sample
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
        self.obs_key_age = x  # (optional, see above, do not provide if .age is provided)
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
        # This cell type annotation is free text but is mapped to an ontology via a .csv file with the same name and
        # directory as the python file of this data loader (see below).


2. A function called to load the data set into memory:
It is important to set an automated path indicating the location of the raw files here.
Our recommendation for this directory set-up is that you define a directory folder in your directory structure
in which all of these raw files will be (self.path) and then add a sub-directory named as
`self.directory_formatted_doi` (ie. the doi with all special characters replaced by "_" and place the raw files
directly into this sub directory.

.. code-block:: python

    def _load(self, fn=None):
        # assuming that i uploaded an h5ad somewhere (in self.download)
        if fn is None:
            fn = os.path.join(self.path, self.directory_formatted_doi, "my.h5ad")
        self.adata = anndata.read(fn)  # loading instruction into .adata, use other ones if the data is not h5ad
        # Some times, you need to load multiple files (e.g. counts and annotation), all of this code would be here.


In summary, a simply example data loader for a mouse lung data set could look like this:

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
            self.doi = "my preprint"
            self.download_url_data = "my GEO upload"
            self.normalisation = "raw"  # because I uploaded raw counts, which is good practice!
            self.organ = "lung"
            self.organism = "mouse"
            self.assay_sc = "smart-seq2"
            self.year = "2020"
            self.sample_source = "primary_tissue"

            self.obs_key_cellontology_original = "louvain_named"  # i save my cell type names in here

        def _load(self, fn=None):
            # assuming that i uploaded an h5ad somewhere (in self.download)
            if fn is None:
                fn = os.path.join(self.path, self.directory_formatted_doi, "my.h5ad")
            self.adata = anndata.read(fn)


Data loaders can be added into a copy of the sfaira repository and can be used locally before they are contributed to
the public sfaira repository.
Alternatively, we also provide the optional dependency sfaira_extensions (https://github.com/theislab/sfaira_extension)
in which local data and cell type annotation can be managed separately but still be loaded as usual through sfaira.
The data loaders and cell type annotation formats between sfaira and sfaira_extensions are identical and can be easily
copied over.

Handling multiple data sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have multiple data sets in a study which are all saved in separate files which come in similar formats:
You can subclass `DatasetBaseGroupLoadingManyFiles` instead of `DatasetBase` and proceed as usual,
only with adding `SAMPLE_FNS` in the data loader file name space,
which is a list of all file names addressed with this file.
You can then refer to an additional property of the Dataset class, `self.sample_fn` during loading
or when dynamically defining meta data in the constructor.
Note that you can always add additional data loaders for further, less streamlined, data sets to such a study.

If you have multiple data sets in a study which are all saved in one file:
You can subclass `DatasetBaseGroupLoadingOneFile` instead of `DatasetBase` and proceed as usual,
only with adding `SAMPLE_IDS` in the data loader file name space,
which is a list of all sample IDs addressed with this file.
You can then refer to an additional property of the Dataset class, `self.sample_id` during loading
or when dynamically defining meta data in the constructor.
Note that `self.sample_id` refers to a `self.adata.obs` column in the loaded data set,
this column has to be defined in `self.obs_key_sample`, which needs to be defined in the constructor.
Note that you can always add additional data loaders for further, less streamlined, data sets to such a study.

Creating dataloaders with the commandline interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sfaira features an interactive way of creating, formatting and testing dataloaders.
The common workflow look as follows:

1. Create a new dataloader with ``sfaira create-dataloader``
2. Validate the dataloader with ``sfaira lint-dataloader <path>``

When creating a dataloader with ``sfaira create-dataloader`` common information such as
your name and email are prompted for, followed by dataloader specific attributes such as organ, organism and many more.
If the requested information is not available simply hit enter and continue until done. If you have mixed organ or organism
data you will have to resolve this manually later. Your dataloader template will be created in your current working directory
in a folder resembling your doi.

The created files are:

.. code-block::

    ├── extra_description.txt <- Optional extra description file
    ├── __init__.py
    ├── NA_NA_2021_NA_Einstein_001.py <- Contains the load function to load the data
    ├── NA_NA_2021_NA_Einstein_001.yaml <- Specifies all data loader data

Now simply fill in all missing properties in your dataloader scripts and yaml file.
When done optionally run ``sfaira clean-dataloader <path to *.yaml>`` on the just filled out dataloader yaml file.
All unused attributes will be removed.

Next validate the integrity of your dataloader content with ``sfaira lint-dataloader <path to *.yaml>``.
All tests must pass! If any of the tests fail please revisit your dataloader and add the missing information.

Map cell type labels to ontology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The entries in `self.obs_key_cellontology_original` are free text but are mapped to an ontology via a .csv file with
the same name and directory as the python file in which the data loader is located.
This .csv contains two columns with one row for each unique cell type label.
The free text identifiers in the first column "source",
and the corresponding ontology term in the second column "target".
You can write this file entirely from scratch.
Sfaira also allows you to generate a first guess of this file using fuzzy string matching
which is automatically executed when you run the template data loader unit test for the first time with you new loader.
Conflicts are not resolved in this first guess and you have to manually decide which free text field corresponds to which
ontology term in the case of conflicts.
Still, this first guess usually drastically speeds up this annotation harmonization.

Cell type ontology management
-----------------------------

Sfaira maintains a wrapper of the Cell Ontology as a class which allows additions to this ontology.
This allows us to use the core ontology used in the community as a backbone and to keep up with newly identifed cell types on our own.
We require all extensions of the core ontology not to break the directed acyclic graph that is the ontology:
Usually, such extensions would be additional leave nodes.

Second, we maintain cell type universes for anatomic structures.
These are dedicated for cell type-dependent models which require a defined set of cell types.
Such a universe is a set of nodes in the ontology.

Contribute cell types to ontology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please open an issue on the sfaira repo with a description what type of cell type you want to add.

Using ontologies to train cell type classifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cell type classifiers can be trained on data sets with different coarsity of cell type annotation using aggregate
cross-entropy as a loss and aggregate accuracy as a metric.
The one-hot encoded cell type label matrix is accordingly modified in the estimator class in data loading if terms
that correspond to intermediate nodes (rather than leave nodes) are encountered in the label set.

Metadata management
-------------------

We constrain meta data by ontologies where possible. The current restrictions are:

    - .age: unconstrained string, try using units of years for human, units of months for mice and units of days for
        cell culture samples
    - .dev_stage: unconstrained string, this will constrained to an ontology in the future,
        try choosing from HSAPDV (http://www.obofoundry.org/ontology/hsapdv.html) for human
        or from MMUSDEV (http://www.obofoundry.org/ontology/mmusdv.html) for mouse
    - .cell_line: unconstrained string, this will be constrained to an ontology later. try choosing from cellosaurus
        cell line database (https://web.expasy.org/cellosaurus/)
    - .ethnicity: unconstrained string, this will constrained to an ontology in the future,
        try choosing from HANCESTRO (https://www.ebi.ac.uk/ols/ontologies/hancestro)
    - .healthy: bool
    - .normalisation: unconstrained string, this will constrained to an ontology in the future,
        try using {"raw", "scaled"}
    - .organ: unconstrained string, this will constrained to an ontology in the future, try to choose
        term from Uberon (http://www.obofoundry.org/ontology/ehdaa2.html)
        or from EHDAA2 (http://www.obofoundry.org/ontology/ehdaa2.html) for human
        or from EMAPA (http://www.obofoundry.org/ontology/emapa.html) for mouse
    - .organism: constrained string, {"mouse", "human"}. In the future, we will use NCBITAXON
        (http://www.obofoundry.org/ontology/ncbitaxon.html).
    - .assay_sc: unconstrained string, this will constrained to an experimental protocol ontology in the future,
        try choosing a term from https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0010183&viewMode=All&siblings=false
    - .assay_differentiation: unconstrained string, try to provide a base differentiation protocol (eg. Lancaster, 2014)
        as well as any amendments to the original protocol
    - .assay_type_differentiation: constrained string, {"guided", "unguided"}
    - .sample_source: constrained string, {"primary_tissue", "2d_culture", "3d_culture", "cancer"}
    - .sex: constrained string, {"female", "male"}
    - .state_exact: unconstrained string, try to be concise and anticipate that this field is queried by automatised searches.
        If you give treatment concentrations, intervals or similar measurements use square brackets around the quantity
        and use units: `[1g]`
    - .year: must be an integer year, e.g. 2020

Follow this issue_ for details on upcoming ontology integrations.

.. _issue: https://github.com/theislab/sfaira/issues/16
