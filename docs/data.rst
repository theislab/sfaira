Data
======

.. image:: https://raw.githubusercontent.com/theislab/sfaira/master/resources/images/data_zoo.png
   :width: 600px
   :align: center

Build data repository locally
------------------------------

Build a repository structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    1. Choose a directory to dedicate to the data base, called root in the following.
    2. Run the sfaira download script (sfaira.data.utils.download_all). Alternatively, you can manually set up a data base by making subfolders for each study.

Note that the automated download is a feature of sfaira but not the core purpose of the package:
Sfaira allows you efficiently interact with such a local data repository.
Some data sets cannot be automatically downloaded and need you manual intervention, which we report in the download script output.

Use 3rd party repositories
~~~~~~~~~~~~~~~~~~~~~~~~~~
Some organization provide streamlined data objects that can be directly consumed by data zoos such as sfaira.
One example for such an organization is the cellxgene_ data portal.
Through these repositories, one can easily build or extend a collection of data sets that can be easily interfaced with sfaira.
Data loaders for cellxgene structured data objects will be available soon!
Contact us for support of any other repositories.

.. _cellxgene: https://cellxgene.cziscience.com/

Add data sets
~~~~~~~~~~~~~

    1. Write a data loader as outlined below.
    2. Identify the raw files as indicated in the data loader classes and copy them into your directory structure as required by your data loader.
    3. You can contribute the data loader to public sfaira, we do not manage data upload though. During publication, you would upload this data set to a server like GEO and the data loader contributed to sfaira would use this download link.

Use data loaders on existing data repository
--------------------------------------------

You only want to use data sets with existing data loaders and have adapted your directory structure as above?
In that case, you can immediately start using the data loader functions, you just need to supply the root directory
of the directory structure as `path to the constructor of the class that you are using.
Depending on the functionalities you want to use, you would often want to create a directory with cached meta data
first. This can be easily done via the script sfaira.data.utils.create_meta.py. This meta information is necessary to
anticipate file sizes for backing merged adata objects, for example, and is used for lazy loading.

Write data loaders
------------------

The study-centric data loader module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the sfaira code, data loaders are organised into directories, which correspond to publications.
All data loaders corresponding to data sets of one study are grouped into this directory.
This directory contains an `__init__.py` file which makes these data loaders visible to sfaira:

.. code-block:: python

    FILE_PATH = __file__


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
preprint and publication DOIs if both are available.
We will also mention publication names in issues, you will however not find these in the code.

.. _code: https://github.com/theislab/sfaira/tree/dev
.. _issues: https://github.com/theislab/sfaira/issues


The data loader python file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each data set (organsism, organ, protocol, optionally also batches) has its own data loader class. Each such class is
in a separate file and inherits from a base class that contains most functionalities. Accordingly, the data loader class
looks very similar in parts to a cell in a juypter notebook that performs data loading. We suggest to copy a data loader
class file and simply adapt to the new data. The core features that must be included are:

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

        self.id = x  # unique identifier of data set (Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi).

        self.author = x  # author (list) who sampled / created the data set
        self.doi = x  # doi of data set accompanying manuscript

        self.download_url_data = x  # download website(s) of data files
        self.download_url_meta = x  # download website(s) of meta data files

        self.age = x  # (*, optional) age of sample
        self.dev_stage = x  # (*, optional) developmental stage of organism
        self.ethnicity = x  # (*, optional) ethnicity of sample
        self.healthy = x  # (*, optional) whether sample represents a healthy organism
        self.normalisation = x  # (optional) normalisation applied to raw data loaded (ideally counts, "raw")
        self.organ = x  # (*, optional) organ (anatomical structure)
        self.organism = x  # (*) species / organism
        self.protocol = x  # (*, optional) protocol used to sample data (e.g. smart-seq2)
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
        self.obs_key_protocol = x  # (optional, see above, do not provide if .protocol is provided)
        self.obs_key_sex = x  # (optional, see above, do not provide if .sex is provided)
        self.obs_key_state_exact = x  # (optional, see above, do not provide if .state_exact is provided)
        # Additionally, cell type annotation is ALWAYS provided per cell in .obs, this annotation is optional though.
        # name of column which contain streamlined cell ontology cell type classes:
        self.obs_key_cellontology_original = x  # (optional)
        # This cell type annotation is free text but is mapped to an ontology via a .csv file with the same name and
        # directory as the python file of this data loader (see below).

        # A dictionary of dictionaries with:
        # One item for each annotation label that is not contained in the ontology.
        # This item maps a custom ID to an ontology supported ID.
        # Note that you have to load your custom IDs, to which this refers to, in load().
        self.class_maps = {
            "0": {  # one entry for each cell type version for this species and organ
                'my weird name for T cells': 'T cell',  # one map from a custom ID to an ontology supported ID
            },
        }


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
            self.protocol = "smart-seq2"
            self.year = "2020"

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

Map cell type labels to ontology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The entries in `self.obs_key_cellontology_original` are free text but are mapped to an ontology via a .csv file with
the same name and directory as the python file in which the data loader is located.
This .csv contains two columns with one row for each unique cell type label and their free text identifiers in the first
column, and the corresponding ontology term in the second column.
You could write this file entirely from scratch.
Sfaira also allows you to generate a first guess of this file using fuzzy string matching via ToDo.
Conflicts are not resolved in this first guess and you have to manually decide which free text field corresponds to which
ontology term in the case of conflicts.
Still, this first guess usually drastically speeds up this annotation harmonization.


Repetitive data loader code
~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are instances in which you find yourself copying code between data loader files corresponding to one study.
In most of these cases, you can avoid the copy operations and share the code more efficiently.

If you have multiple data files which each correspond to a data set and are structured similarly, you can define a super
class which contains the shared constructor and `_load()` code, from which each data set specific loader inherits.
ToDo: Example.

If you have a single file which contains the data from multiple data sets which belong to a data loader each,
because of different meta data or batches for example,
you can set up a `group.py` file which defines a DatasetGroup for this study, which controls the generation of Datasets.
ToDo: Example.

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

    - .age: unconstrained string, try using units of years for human and units of months for mice
    - .dev_stage: unconstrained string, this will constrained to an ontology in the future,
        try choosing from HSAPDV (http://www.obofoundry.org/ontology/hsapdv.html) for human
        or from MMUSDEV (http://www.obofoundry.org/ontology/mmusdv.html) for mouse
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
    - .protocol: unconstrained string, this will constrained to an anatomic ontology in the future,
        try choosing a term from https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0010183&viewMode=All&siblings=false
    - .sex: constrained string, {"female", "male"}
    - .state_exact: unconstrained string, try to be concise and anticipate that this field is queried by automatised searches.
        If you give treatment concentrations, intervals or similar measurements use square brackets around the quantity
        and use units: `[1g]`
    - .year: must be an integer year, e.g. 2020

Follow this issue_ for details on upcoming ontology integrations.

.. _issue: https://github.com/theislab/sfaira/issues/16

Genome management
-----------------

You do not have to worry about this unless you are interested,
this section is not required reading for writing data loaders.

We streamline feature spaces used by models by defining standardized gene sets that are used as model input.
Per default, sfaira works with the protein coding genes of a genome assembly right now.
A model topology version includes the genome it was trained for, which also defines the feature of this model as genes.
As genome assemblies are updated, model topology version can be updated and models retrained to reflect these changes.
Note that because protein coding genes do not change drastically between genome assemblies,
sample can be carried over to assemblies they were not aligned against by matching gene identifiers.
Sfaira automatically tries to overlap gene identifiers to the genome assembly selected through the current model.

FAQ
---

How is the dataset’s ID structured?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi

How do I assemble the data set ID if some of its element meta data are not unique?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The data set ID is designed to be a unique identifier of a data set.
Therefore, it is not an issue if it does not capture the full complexity of the data.
Simply choose the meta data value out of the list of corresponding values which comes first in the alphabet.

What are cell-wise and sample-wise meta data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Metadata can be set on a per sample level or, in some cases, per cell.
Sample-wise meta data can be directly set in the constructor (e.g self.organism = “human”).
Cell-wise metadata can be provided in `.obs` of the loaded data, here,
a Dataset attribute contains the name of the `.obs` column that contains these cell-wise labels
(e.g. self.obs_key_organism).
Note that sample-wise meta data should be yielded as such and not as a column in `.obs` to simplify loading.

Which meta data objects are optional?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Mandatory on sample (self.attribute) or cell level (self.obs_key_attribute):

    - .id: Dataset ID. This is used to identify the data set uniquely.
        Example: self.id = "human_colon_2019_10x_smilie_001_10.1016/j.cell.2019.06.029"
    - .download_url_data: Link to data download website.
        Example: self.download = "some URL"
    - .download_url_meta: Download link to metadata. Assumes that meta data is defined in .download_url_data if not
        specified.
        Example: self.download_meta = "some URL"
    - .var_symbol_col, .var_ensembl_col: Location of gene name as gene symbol and/or ENSEMBL ID in adata.var
        (if index of adata.var, set to “index”, otherwise to column name). One of the two must be provided.
        Example: self.var_symbol_col = 'index', self.var_ensembl_col = “GeneID”
    - .author: First author of publication (or list of all authors).
        self.author = "Last name, first name" # or ["Last name, first name", "Last name, first name"]
    - .doi: Doi of publication
        Example: self.doi = "10.1016/j.cell.2019.06.029"
    - .organism (or .obs_key_organism): Organism sampled.
        Example: self.organism = “human”

Highly recommended:

    - .normalization: Normalization of count data:
        Example: self.normalization = “raw”
    - .organ (or .obs_key_organ): Organ sampled.
        Example: self.organ = “liver”
    - .protocol (or .obs_key_protocol): Protocol with which data was collected.
        Example: self.protocol = “10x”

Optional (if available):

    - .age (or .obs_key_age): Age of individual sampled.
        Example: self.age = 80  # (80 years old for human)
    - .dev_stage (or .obs_key_dev_stage): Developmental stage of individual sampled.
        Example: self.dev_stage = “mature”
    - .ethnicity (or .obs_key_ethnicity): Ethnicity of individual sampled (only for human).
        Example: self.ethnicity = “free text”
    - .healthy (or .obs_key_healthy): Is the sampled from a disease individual? (bool)
        Example: self.healthy = True
    - .sex (or .obs_key_sex): Sex of individual sampled.
        Example: self.sex = “male”
    - .state_exact (or .obs_key_state_exact): Exact disease state
        self.state_exact = free text
    - .obs_key_cellontology_original: Column in .obs in which free text cell type names are stored.
        Example: self.obs_key_cellontology_original = 'CellType'
    - .year: Year of publication:
        Example: self.year = 2019

How do I cache data sets?
~~~~~~~~~~~~~~~~~~~~~~~~~
When loading a dataset with `Dataset.load(),`you can specify if the adata object
should be cached or not  (allow_caching= True).
If set to True, the loaded adata object will be cached as an h5ad object for faster reloading.

How do I add cell type annotation?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We are simplifying this right now, new instructions will be available second half of January.

Why are constructor (`__init__`) and loading function (`_load`) split in the template data loader?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Initiation and data set loading are handled separately to allow lazy loading.
All steps that are required to load the count data and
additional metadata should be defined solely in the `_load` section.
Setting of class metadata such as `.doi`, `.id` etc. should be done in the constructor.

How do I tell sfaira where the gene names are?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
By setting the attributes `.var_symbol_col` or `.var_ensembl_col` in the constructor.
If the gene names are in the index of this data frame, you can set “index” as the value of these attributes.

I only have gene symbols (human readable names, often abbreviations), such as HGNC or MGI, but not ENSEMBL identifiers, is that a problem?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No, that is not a problem. They will automatically be converted to Ensembl IDs.
You can, however, specify the reference genome in `Dataset.load(match_to_reference = ReferenceGenomeName)`
to which the names should be mapped to.

I have CITE-seq data, where can I put the protein quantification?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will soon provide a structured interface for loading and accessing CITE-seq data,
for now you can add it into `self.adata.obsm[“CITE”]`.
