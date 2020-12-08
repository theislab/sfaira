.. module:: sfaira
.. automodule:: sfaira
   :noindex:

API
===

Import sfaira as::

   import sfaira.api as sfaira



Data: `data`
------------

.. module:: sfaira.data
.. currentmodule:: sfaira

The sfaira data zoo API.


Pre-defined data set collections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This sub-module gives you access to curated subsets of the data zoo, e.g. all data sets from human lungs.

.. autosummary::
   :toctree: .

   data.human
   data.mouse


Functionalities for interactive data analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This sub-module gives you access to functionalities you need to define your own data set collections based on the sfaira data zoo.

.. autosummary::
   :toctree: .

   data.DatasetBase
   data.DatasetGroupBase
   data.DatasetSuperGroup


Functionalities for interactive data analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This sub-module gives you access to functionalities you need to load new data live into the data zoo to handle a raw data set in the context of zoo data sets.

.. autosummary::
   :toctree: .

   data.DatasetInteractive


Genomes: `genomes`
------------------

.. module:: sfaira.genomes
.. currentmodule:: sfaira

This sub-module gives you access to properties of the genome representations used in sfaira.

.. autosummary::
   :toctree: .

   genomes.ExtractFeatureListEnsemble


Models: `models`
----------------

.. module:: sfaira.models
.. currentmodule:: sfaira

The sfaira model zoo API for advanced use.
This API is structured by streamlined, task-specific APIs for specific analysis problems.
This API is targeted at developers, see also `ui` for a user centric wrapping API for this model zoo.


Cell-type predictor models
~~~~~~~~~~~~~~~~~~~~~~~~~~

This sub-module handles models that predict cell types.

.. autosummary::
   :toctree: .

   models.celltype


Embedding models
~~~~~~~~~~~~~~~~

This sub-module handles models that embed expression vectors (cells) into a latent space.

.. autosummary::
   :toctree: .

   models.embedding


Train: `train`
--------------

.. module:: sfaira.train
.. currentmodule:: sfaira

The interface for training sfaira compatible models.
This is a sub-module dedicated for developers to ease model training and deployment.

Trainer classes
~~~~~~~~~~~~~~~

Trainer class wrap estimator classes (which wrap model classes) and handle grid-search specific tasks centred on model fits,
such as saving evaluation metrics and model weights.

.. autosummary::
   :toctree: .

   train.TargetZoos
   train.TrainModelCelltype
   train.TrainModelEmbedding


Grid search summary classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grid search summary classes allow a developer to easily interact with a finished grid search by loading and summarising results,
which were saved through Trainer classes.

.. autosummary::
   :toctree: .

   train.GridsearchContainer
   train.SummarizeGridsearchCelltype
   train.SummarizeGridsearchEmbedding

User interface: `ui`
--------------------

.. module:: sfaira.ui
.. currentmodule:: sfaira

This sub-module gives users access to the model zoo, including model query from remote servers.
This API is designed to be used in analysis workflows and does not require any understanding of the way models are defined and stored.

.. autosummary::
   :toctree: .

   ui.UserInterface
