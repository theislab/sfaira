.. module:: sfaira
.. automodule:: sfaira
   :noindex:

API
====

Import sfaira as::

   import sfaira



Data: `data`
-------------

.. module:: sfaira.data
.. currentmodule:: sfaira

The sfaira data zoo API.

Dataset representing classes used for development:

.. autosummary::
   :toctree: .

   data.DatasetBase
   data.DatasetGroup
   data.DatasetGroupDirectoryOriented
   data.DatasetSuperGroup

Interactive data class to use a loaded data object in the context sfaira tools:

.. autosummary::
   :toctree: .

   data.DatasetInteractive

Dataset universe to interact with all data loader classes:

.. autosummary::
   :toctree: .

   data.Universe

Data store handling:

.. autosummary::
   :toctree: .

   data.load_store
   data.DistributedStoreBase
   data.DistributedStoreDao
   data.DistributedStoreH5ad


Estimator classes: `estimators`
--------------------------------

.. module:: sfaira.estimators
.. currentmodule:: sfaira

Estimator classes from the sfaira model zoo API for advanced use.

.. autosummary::
   :toctree: .

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
   :toctree: .

   models.celltype.CellTypeMarker
   models.celltype.CellTypeMarker
   models.celltype.CellTypeMlp
   models.celltype.CellTypeMlpVersioned

Embedding models
~~~~~~~~~~~~~~~~~
Classes that wrap tensorflow embedding models.

.. autosummary::
   :toctree: .

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
   :toctree: .

   train.TrainModelCelltype
   train.TrainModelEmbedding

Grid search summaries
~~~~~~~~~~~~~~~~~~~~~
Classes to pool evaluation metrics across fits in a grid search.

.. autosummary::
   :toctree: .

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
   :toctree: .

   versions.genomes.GenomeContainer

Metadata
~~~~~~~~~
Dataset metadata management.
Base classes to manage ontology files:

.. autosummary::
   :toctree: .

   versions.metadata.Ontology
   versions.metadata.OntologyList
   versions.metadata.OntologyHierarchical
   versions.metadata.OntologyObo
   versions.metadata.OntologyOboCustom

Onotology-specific classes:

.. autosummary::
   :toctree: .

   versions.metadata.OntologyCellosaurus
   versions.metadata.OntologyCl
   versions.metadata.OntologyHsapdv
   versions.metadata.OntologyMondo
   versions.metadata.OntologyMmusdv
   versions.metadata.OntologySinglecellLibraryConstruction
   versions.metadata.OntologyUberon

Class wrapping cell type ontology for predictor models:

.. autosummary::
   :toctree: .

   versions.metadata.CelltypeUniverse

Topologies
~~~~~~~~~~~
Model topology management.

.. autosummary::
   :toctree: .

   versions.topologies.TopologyContainer

User interface: `ui`
---------------------

.. module:: sfaira.ui
.. currentmodule:: sfaira

This sub-module gives users access to the model zoo, including model query from remote servers.
This API is designed to be used in analysis workflows and does not require any understanding of the way models are defined and stored.

.. autosummary::
   :toctree: .

   ui.UserInterface
