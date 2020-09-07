Ontologies
==========

Introduction to sfaira ontology management
------------------------------------------

Sfaira maintains versioned cell type universes and ontologies by species and organ.
A cell type universe is a list of the unique, most fine-grained cell type definitions available.
These cell types can be referred to by a human readable cell type name or a structure identifier within an ontology,
an ontology ID.
Often, one is also interested in access to more coarse grained groups of cell types, for example if the data quality
does not allow to distinguish between T cell subtypes.
To allow coarser type definition, sfaira maintains hierarchies of cell types, in which each hierarchical level is again
defined by a cell type identifier.
Such a hierarchy can be writted as directed acyclic graph which has the cell type universe as its leave nodes.
Intuitively, the cell type hierarchy graph depends on the cell type universe.
Accordingly, both are versioned together in sfaira:
Updates in the cell type universe, such as discovery of a new cell type, lead to an update of the ontology and an
incrementation in both of their versions.
These versioned changes materialise as a distinct list (universe) and dictionary (ontology) for each version in the
file that harbors the species- and organ-specific class that inherits from CelltypeVersionsBase and thus are available
even after updates.
This versioning without depreceation of the old objects allows sfaira to execute and train models that were designed
for older cell type universes and thus ensures reproducibility.

Contribute cell types to ontologies
-----------------------------------

To contibute new cell types or change existing cell type universe entries, the cell type universe version has to be
incremented and the new entry can simply be added to the list or modified in the list.
We do not increment the universe version if a change does not influence the identity of a leave node with respect to
the other types in the universe, ie if it simply changes the spelling of a cell type or if an onology ID is added to
a type that previously did not have one.

Contribute hierarchies to ontologies
------------------------------------

To contribute a term to a cell type ontology, one just has to add a dictionary item that defines the new term as a set
of the leave nodes (cell type universe) of the corresponding universe version.


Using ontologies to train cell type classifiers
-----------------------------------------------

Cell type classifiers can be trained on data sets with different coarsity of cell type annotation using aggregate
cross-entropy as a loss and aggregate accuracy as a metric.
The one-hot encoded cell type label matrix is accordingly modified in the estimator class in data loading if terms
that correspond to intermediate nodes (rather than leave nodes) are encountered in the label set.
