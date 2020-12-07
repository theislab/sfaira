Environment
===========

scanpy
------

scanpy_ provides an environment of tools that can be used to analysis single-cell data in python.
sfaira allows users to easily query third party data sets and models to complement these analysis workflows.

.. _scanpy: https://github.com/theislab/scanpy

Study-centric data set servers
------------------------------

Many authors of data sets provide their data sets on servers such as GEO or in manuscript supplements.
Our data zoo interface is able to represent these data sets such that they can be queried in a streamlined fashion,
together with many other data sets.

Data providers which streamline data
------------------------------------

Some organization provide streamlined data objects that can be directly consumed by data zoos such as sfaira.
One example for such an organization is the cellxgene_ data portal.
Through these repositories, one can easily build or extend a collection of data sets that can be easily interfaced with sfaira.
Data loaders for cellxgene structured data objects will be available soon!
Contact us for support of any other repositories.

.. _cellxgene: https://cellxgene.cziscience.com/


Single-cell study data base by Svensson et al.
----------------------------------------------

Svensson_ et al. published a single-cell database_ in the form of a table in which each row contains a description of a study which published single-cell RNA-seq data.
Some of these data sets are already included in sfaira:
We aim to include all of these data sets but have currently a much smaller set of data sets.
However, the sfaira database handles direct access to these data sets,
allowing users to interact with anndata_ objects that represent subsets of this database,
without requiring any data loading code from the user,
therefore opening up further use cases.
Consider also our interactive website_ for a graphical user interface to our complete data zoo.

.. _Svensson: https://academic.oup.com/database/article/doi/10.1093/database/baaa073/6008692
.. _database: https://www.nxn.se/single-cell-studies/gui
.. _anndata: https://github.com/theislab/anndata
.. _website: https://theislab.github.io/sfaira-site/index.html
