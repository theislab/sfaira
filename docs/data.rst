Data
======

Build data repository locally
------------------------------

Build a repository structure:
1. Choose a directory to dedicate to the data base, called root in the following.
2. Make subfolders in root for each organism for which you want to build a data base.
3. Make subfolders for each organ whithin each organism for which you want to build a data base.

Add data sets:
4. For each species and organ combination, choose the data sets that you want to use.
5. Identify the raw files as indicated in the data loader classes and copy them into the folder. Use processed data
using the described processing if this is required: This is usually done to speed up loading for file
formats that are difficult to access.

Use data loaders on existing data repository
--------------------------------------------

You only want to use data sets with existing data loaders and have adapted your directory structure as above?
In that case, you can immediately start using the data loader functions, you just need to supply the root directory
of the directory structure as `path to the constructor of the class that you are using.
Depending on the functionalities you want to use, you need to create a directory with data set meta data first. This
can be easily done via the data set api itself, example python scripts are under benchmarks/data_preparation. This
meta information is necessary to anticipate file sizes for backing merged adata objects for example.

TODO example.

Contribute data loaders
-----------------------

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
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = x  # your-species
        self.id = x  # "organism_organ_year_protocoll_first-author_doi"
        self.download_website = x  # link to raw data
        self.organ = x  #y ourorgan
        self.sub_tissue = x # sub-tissue name, otherwise organ
        self.dev_stage = x  # developmental stage of organism
        self.has_celltypes = x  # if cell type annotation is available

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
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "human/eye/my_data.h5ad")  defined file in streamlined directory structure
        self.adata = anndata.read(fn)  # loading instruction into .adata, use other ones if the data is not h5ad

        self.adata.uns["lab"] = x  # load the adata.uns with meta data
        self.adata.uns["year"] = x
        self.adata.uns["doi"] = x
        self.adata.uns["protocol"] = x  # e.g. 10x, microwell, seqwell...
        self.adata.uns["organ"] = self.organ
        self.adata.uns["subtissue"] = self.sub_tissue
        self.adata.uns["animal"] = x
        self.adata.uns["id"] = self.id
        self.adata.uns["wget_download"] = self.download_website
        self.adata.uns["has_celltypes"] = self.has_celltypes
        self.adata.uns["counts"] = 'raw'
        self.adata.uns["dev_stage"] = self.dev_stage

        # Class expects unprocessed cell type labels in self.adata.obs["cell_ontology_class"]
        self.adata.obs["cell_ontology_class"] = self.adata.obs['CellType']
        # You can additional set self.adata.obs["cell_ontology_id"] if you have streamlined ontology IDs. This are also
        # defined in the cell type universe lists.
        self.adata.obs["healthy"] = x  # boolean tissue sample healthy or diseased / treated
        self.adata.obs["state_exact"] = x  # exact tissue state as string, e.g. "tumor" or "healthy"

        self._convert_and_set_var_names(symbol_col='names', ensembl_col='ensembl', new_index='ensembl')



Data loaders can be added into a copy of the sfaira repository and can be used locally before they are contributed to
the public sfaira repository.
Alternatively, we also provide the optional dependency sfaira_extensions (https://github.com/theislab/sfaira_extension)
in which local data and cell type annotation can be managed separately but still be loaded as usual through sfaira.
The data loaders and cell type annotation formats between sfaira and sfaira_extensions are identical and can be easily
copied over.


Handling ontologies in data loaders
-----------------------------------

Each data loader has a versioned cell type annotation map, a dictionary.
This dictionary allows mapping of the cell type annotations that come with the raw form of the data set to the cell type
universe or ontology terms defined in sfaira, this is, however, only done upon loading of the data (.load()).
The outcome of this map is a new set of cell type labels that can be propagated to leave nodes of the ontology graph.
This dictionary requires a new entry for each new version of the corresponding cell type universe.

