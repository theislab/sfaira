FAQ
====

Data zoo
---------

How is load() function used in data loading?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`load()` contains all processing steps that load raw data files into a ready to use adata object.
This adata object can be cached as an h5ad file named after the dataset ID for faster reloading
(if allow_caching=True), which exactly skips the code in `load()`.
`load()` can be triggered to reload from scratch even if cached data is available
(if use_cached=False).

How is the feature space (gene names) manipulated during data loading?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sfaira provides both gene names and ENSEMBL IDs. Missing IDs will automatically be inferred from the gene names and
vice versa.
Version tags on ENSEMBL gene IDs will be removed if specified (if remove_gene_version=True);
in this case, counts are aggregated across these features.
Sfaira makes sure that gene IDs in a dataset match IDs of chosen reference genomes.

Datasets, DatasetGroups, DatasetSuperGroups - what are they?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dataset: Custom class that loads a specific dataset.
DatasetGroup: A dataset group manages collection of data loaders (multiple instances of Dataset).
This is useful to group for example all data loaders corresponding to a certain study or a certain tissue.
DatasetSuperGroups: A group of DatasetGroups that allow easy addition of multiple instances of DatasetGroup.

Basics of sfaira lazy loading via split into constructor and load() function.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The constructor of a dataset defines all metadata associated with this data set.
The loading of the actual data happens in the `load()` function and not in the constructor.
This is useful as it allows initialising the datasets and accessing dataset metadata
without loading the actual count data.
DatasetGroups can contain initialised Datasets and can be sub-setted based on metadata
before loading is triggered across the entire group.
