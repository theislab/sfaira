API
====

Import sfaira as::

   import sfaira



Data: `data`
-------------

Data loaders
~~~~~~~~~~~~~

.. module:: sfaira.data
.. currentmodule:: sfaira

The sfaira data zoo API.

Dataset representing classes used for development:

.. autosummary::
   :toctree: api

   data.DatasetBase
   data.DatasetGroup
   data.DatasetGroupDirectoryOriented
   data.DatasetSuperGroup

Interactive data class to use a loaded data object in the context sfaira tools:

.. autosummary::
   :toctree: api

   data.DatasetInteractive

Dataset universe to interact with all data loader classes:

.. autosummary::
   :toctree: api

   data.Universe

Stores
~~~~~~~

.. module:: sfaira.data
.. currentmodule:: sfaira

We distinguish stores for a single feature space, which could for example be a single organism,
and those for multiple feature spaces.
Critically, data from multiple feature spaces can be represented as a data array for each feature space.
In `load_store` we represent a directory of datasets as a instance of a multi-feature space store and discover all feature
spaces present.
This store can be subsetted to a single store if only data corresponding to a single organism is desired,
for example.
The core API exposed to users is:

.. autosummary::
   :toctree: api

   data.load_store

Store classes for a single feature space:

.. autosummary::
   :toctree: api

   data.DistributedStoreSingleFeatureSpace

Store classes for a multiple feature spaces:

.. autosummary::
   :toctree: api

   data.DistributedStoreMultipleFeatureSpaceBase
   data.DistributedStoresAnndata
   data.DistributedStoresDao
   data.DistributedStoresH5ad

Carts
~~~~~~

.. module:: sfaira.data
.. currentmodule:: sfaira

Stores represent on-disk data collection and perform operations such as subsetting.
Ultimatively, they are often used to emit data objects, which are "carts".
Carts are specific to the underlying store's data format and expose iterators, data matrices and
adaptors to machine learning framework data pipelines, such as tensorflow and torchc data.
Again, carts can cover one or multiple feature spaces.

.. autosummary::
   :toctree: api

   data.store.carts.CartSingle
   data.store.carts.CartMulti

The emission of data from cart iterators and adaptors is controlled by batch schedules,
which direct how data is released from the underlying data matrix:

.. autosummary::
   :toctree: api

   data.store.batch_schedule.BatchDesignBase
   data.store.batch_schedule.BatchDesignBasic
   data.store.batch_schedule.BatchDesignBalanced
   data.store.batch_schedule.BatchDesignBlocks
   data.store.batch_schedule.BatchDesignFull

For most purposes related to stochastic optimisation, `BatchDesignBasic` is chosen.

Estimator classes: `estimators`
--------------------------------

.. module:: sfaira.estimators
.. currentmodule:: sfaira

Estimator classes from the sfaira model zoo API for advanced use.

.. autosummary::
   :toctree: api

   estimators.EstimatorKeras
   estimators.EstimatorKerasCelltype
   estimators.EstimatorKerasEmbedding

Model classes: `models`
------------------------

.. module:: sfaira.models
.. currentmodule:: sfaira

Model classes from the sfaira model zoo API for advanced use.

Cell type models
~~~~~~~~~~~~~~~~~
Classes that wrap tensorflow cell type predictor models.

.. autosummary::
   :toctree: api

   models.celltype.CellTypeMarker
   models.celltype.CellTypeMarker
   models.celltype.CellTypeMlp
   models.celltype.CellTypeMlpVersioned

Embedding models
~~~~~~~~~~~~~~~~~
Classes that wrap tensorflow embedding models.

.. autosummary::
   :toctree: api

   models.embedding.ModelKerasAe
   models.embedding.ModelAeVersioned
   models.embedding.ModelKerasVae
   models.embedding.ModelVaeVersioned
   models.embedding.ModelKerasLinear
   models.embedding.ModelLinearVersioned
   models.embedding.ModelKerasVaeIAF
   models.embedding.ModelVaeIAFVersioned
   models.embedding.ModelKerasVaeVamp
   models.embedding.ModelVaeVampVersioned

Train: `train`
---------------

.. module:: sfaira.train
.. currentmodule:: sfaira

The interface for training sfaira compatible models.

Trainer classes
~~~~~~~~~~~~~~~~
Classes that wrap estimator classes to use in grid search training.

.. autosummary::
   :toctree: api

   train.TrainModelCelltype
   train.TrainModelEmbedding

Grid search summaries
~~~~~~~~~~~~~~~~~~~~~
Classes to pool evaluation metrics across fits in a grid search.

.. autosummary::
   :toctree: api

   train.GridsearchContainer
   train.SummarizeGridsearchCelltype
   train.SummarizeGridsearchEmbedding

Versions: `versions`
---------------------

.. module:: sfaira.versions
.. currentmodule:: sfaira

The interface for sfaira metadata management.

Genomes
~~~~~~~~
Genome management.

.. autosummary::
   :toctree: api

   versions.genomes.GenomeContainer

Metadata
~~~~~~~~~
Dataset metadata management.
Base classes to manage ontology files:

.. autosummary::
   :toctree: api

   versions.metadata.Ontology
   versions.metadata.OntologyList
   versions.metadata.OntologyHierarchical
   versions.metadata.OntologyObo
   versions.metadata.OntologyOboCustom

Onotology-specific classes:

.. autosummary::
   :toctree: api

   versions.metadata.OntologyCellosaurus
   versions.metadata.OntologyCl
   versions.metadata.OntologyHsapdv
   versions.metadata.OntologyMondo
   versions.metadata.OntologyMmusdv
   versions.metadata.OntologySinglecellLibraryConstruction
   versions.metadata.OntologyUberon

Class wrapping cell type ontology for predictor models:

.. autosummary::
   :toctree: api

   versions.metadata.CelltypeUniverse

Topologies
~~~~~~~~~~~
Model topology management.

.. autosummary::
   :toctree: api

   versions.topologies.TopologyContainer

User interface: `ui`
---------------------

.. module:: sfaira.ui
.. currentmodule:: sfaira

This sub-module gives users access to the model zoo, including model query from remote servers.
This API is designed to be used in analysis workflows and does not require any understanding of the way models are defined and stored.

.. autosummary::
   :toctree: api

   ui.UserInterface
