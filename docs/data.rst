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
    2. Identify the raw files as indicated in the data loader classes and copy them into your directory structure as required by your data laoder.
    3. You can contribute the data loader to public sfaira, we do not manage data upload though. During publication, you would upload this data set to a server like GEO and the dataloader contributed to sfaira would use this download link.

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
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, **kwargs)
        # Data set meta data: You do not have to include all of these and can simply skip lines corresponding
        # to attritbutes that you do not have access to. These are meta data on a sample level.
        # The meta data attributes labeled with (*) may als be supplied per cell, see below,
        # in this case, if you supply a .obs_key* attribute, you ccan leave out the sample-wise attribute.

        self.id = x  # unique identifier of data set (Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi).

        self.author = x  # author (list) who sampled / created the data set
        self.doi = x  # doi of data set accompanying manuscript

        self.download = x  # download website(s) of data files
        self.download_meta = x  # download website(s) of meta data files

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

.. code-block:: python

    def _load(self, fn=None):
        self.adata = anndata.read(fn)  # loading instruction into .adata, use other ones if the data is not h5ad
        # Some times, you need to load multiple files (e.g. counts and annotation), all of this code would be here.


In summary, a simply example data loader for a mouse lung data set could look like this:

.. code-block:: python

    class MyDataset(DatasetBase)
        def __init__(
                self,
                path: Union[str, None] = None,
                meta_path: Union[str, None] = None,
                **kwargs
        ):
            super().__init__(path=path, meta_path=meta_path, **kwargs)
            self.author = "me"
            self.doi = "my preprint"
            self.download = "my GEO upload"
            self.normalisation = "raw"  # because I uploaded raw counts, which is good practice!
            self.organ = "lung"
            self.organism = "mouse"
            self.protocol = "smart-seq2"
            self.year = "2020"

            self.obs_key_cellontology_original = "louvain_named"  # i save my cell type names in here

        def _load(self, fn=None):
            # assuming that i uploaded an h5ad somewhere (in self.download)
            if fn is None:
                fn = os.path.join(self.path, "mouse", "lung", "my.h5ad")
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

    - .organism must either mouse or human.

Follow this issue_ for details on upcoming ontology integrations.

.. _issue: https://github.com/theislab/sfaira/issues/16

Genome management
-----------------

We streamline feature spaces used by models by defining standardized gene sets that are used as model input.
Per default, sfaira works with the protein coding genes of a genome assembly right now.
A model topology version includes the genome it was trained for, which also defines the feature of this model as genes.
As genome assemblies are updated, model topology version can be updated and models retrained to reflect these changes.
Note that because protein coding genes do not change drastically between genome assemblies,
sample can be carried over to assemblies they were not aligned against by matching gene identifiers.
Sfaira automatically tries to overlap gene identifiers to the genome assembly selected through the current model.
