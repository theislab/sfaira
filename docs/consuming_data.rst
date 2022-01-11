.. _consuming_data_rst:

Data loaderes
===============

.. image:: https://raw.githubusercontent.com/theislab/sfaira/master/resources/images/data_zoo.png
   :width: 600px
   :align: center

For a high-level overview of data management in sfaira, read :ref:`data_life_cycle_rst` first.

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

Genome management
-----------------

We streamline feature spaces used by models by defining standardized gene sets that are used as model input.
Per default, sfaira works with the protein coding genes of a genome assembly right now.
A model topology version includes the genome it was trained for, which also defines the feature of this model as genes.
As genome assemblies are updated, model topology version can be updated and models retrained to reflect these changes.
Note that because protein coding genes do not change drastically between genome assemblies,
sample can be carried over to assemblies they were not aligned against by matching gene identifiers.
Sfaira automatically tries to overlap gene identifiers to the genome assembly selected through the current model.
