Environment
===========

scanpy
------

scanpy_ provides an environment of tools that can be used to analysis single-cell data in python.
sfaira allows users to easily query third party data sets and models to complement these analysis workflows.

.. _scanpy: https://github.com/theislab/scanpy

Data zoo
--------

Data providers which streamline data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some organization provide streamlined data objects that can be directly consumed by data zoos such as sfaira.
Examples for such providers are:
- Human Cell Atlas data portal (HCA DCP_)
- cellxgene_ data portal
- Broad_ institute single cell data portal
- EBI_ single cell expression atlas

Through these repositories, one can easily build or extend a collection of data sets that can be easily interfaced with sfaira.
Data loaders for cellxgene structured data objects will be available soon, we are working on interfacing more such organisations!
Contact us for support of any other repositories.

.. _DCP: https://data.humancellatlas.org/explore/
.. _cellxgene: https://cellxgene.cziscience.com/
.. _Broad: https://singlecell.broadinstitute.org/single_cell
.. _EBI: https://www.ebi.ac.uk/gxa/sc/home


Study-centric data set servers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many authors of data sets provide their data sets on servers such as:
- GEO_
- cloud storage servers
- manuscript supplements

Our data zoo interface is able to represent these data sets such that they can be queried in a streamlined fashion,
together with many other data sets.

.. _GEO: https://www.ncbi.nlm.nih.gov/geo/


Single-cell study look-up tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Svensson_ et al. published a single-cell database_ in the form of a table in which each row contains a description of a study which published single-cell RNA-seq data.
Some of these data sets are already included in sfaira,
consider also our interactive website_ for a graphical user interface to our complete data zoo.
Note that this website can be used as a look-up table but sfaira also allows you to directly load and interact with these data sets.

.. _Svensson: https://academic.oup.com/database/article/doi/10.1093/database/baaa073/6008692
.. _database: https://www.nxn.se/single-cell-studies/gui
.. _website: https://theislab.github.io/sfaira-site/index.html
