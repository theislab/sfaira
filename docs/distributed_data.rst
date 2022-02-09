.. _distributed_data_rst:

Distributed data
================

For a high-level overview of data management in sfaira, read :ref:`data_life_cycle_rst` first.
Sfaira supports usage of distributed data for model training and execution.
The tools are summarized under `sfaira.data.store`.
In contrast to using an instance of AnnData in memory, these tools can be used to use data sets that are saved
in different files (because they come from different studies) flexibly and out-of-core,
which means without loading them into memory.
A general use case is the training of a model on a large set of data sets, subsetted by particular cell-wise meta
data, without creating a merged AnnData instance in memory first.

Build a distributed data repository
-----------------------------------

You can use the sfaira dataset API to write streamlined groups of adata instances to a particular disk locaiton that
then is the store directory.
Some of the array backends used for loading stores can read arrays from cloud servers, such as dask.
Therefore, these store directories can also be on cloud servers in some cases.

Reading from a distributed data repository
------------------------------------------

The core use-case is the consumption of data in batches from a python iterator (a "generator").
In contrast to using the full data matrix, this allows for workflows that never require the full data matrix in memory.
This generators can for example directly be used in tensorflow or pytorch stochastic mini-batch learning pipelines.
The core interface is `sfaira.data.load_store()` which can be used to initialise a store instance that exposes a
generator, for example.
An important concept in store reading is that the data sets are already streamlined on disk, which means that they have
the same feature space for example.

Distributed access optimised (DAO) store
----------------------------------------

The DAO store format is a on-disk representation of single-cell data which is optimised for generator-based access and
distributed access.
In brief, DAO stores optimize memory consumption and data batch access speed.
Right now, we are using zarr and parquet, this may change in the future, we will continue to work on this format using
the project name "dao".
Note that data sets represented as DAO on disk can still be read into AnnData instances in memory if you wish!
