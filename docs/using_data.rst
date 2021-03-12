Using Data
==========

.. image:: https://raw.githubusercontent.com/theislab/sfaira/master/resources/images/data_zoo.png
   :width: 600px
   :align: center

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

FAQ
---

How is the dataset’s ID structured?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Organism_Organ_Year_AssaySc_NumberOfDataset_FirstAuthorLastname_doi

How do I assemble the data set ID if some of its element meta data are not unique?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The data set ID is designed to be a unique identifier of a data set.
Therefore, it is not an issue if it does not capture the full complexity of the data.
Simply choose the meta data value out of the list of corresponding values which comes first in the alphabet.

What are cell-wise and sample-wise meta data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Metadata can be set on a per sample level or, in some cases, per cell.
Sample-wise meta data can be directly set in the constructor (e.g self.organism = “human”).
Cell-wise metadata can be provided in `.obs` of the loaded data, here,
a Dataset attribute contains the name of the `.obs` column that contains these cell-wise labels
(e.g. self.obs_key_organism).
Note that sample-wise meta data should be yielded as such and not as a column in `.obs` to simplify loading.

Which meta data objects are mandatory?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Mandatory on sample (self.attribute) or cell level (self.obs_key_attribute):

    - .id: Dataset ID. This is used to identify the data set uniquely.
        Example: self.id = "human_colon_2019_10x_smilie_001_10.1016/j.cell.2019.06.029"
    - .download_url_data: Link to data download website.
        Example: self.download = "some URL"
    - .download_url_meta: Download link to metadata. Assumes that meta data is defined in .download_url_data if not
        specified.
        Example: self.download_meta = "some URL"
    - .var_symbol_col, .var_ensembl_col: Location of gene name as gene symbol and/or ENSEMBL ID in adata.var
        (if index of adata.var, set to “index”, otherwise to column name). One of the two must be provided.
        Example: self.var_symbol_col = 'index', self.var_ensembl_col = “GeneID”
    - .author: First author of publication (or list of all authors).
        self.author = "Last name, first name" # or ["Last name, first name", "Last name, first name"]
    - .doi: Doi of publication
        Example: self.doi = "10.1016/j.cell.2019.06.029"
    - .organism (or .obs_key_organism): Organism sampled.
        Example: self.organism = “human”
    - .sample_source (or .obs_key_sample_source): Whether data was obtained from primary tissue or cell culture
        Example: self.sample_source = "primary_tissue"

Highly recommended:

    - .normalization: Normalization of count data:
        Example: self.normalization = “raw”
    - .organ (or .obs_key_organ): Organ sampled.
        Example: self.organ = “liver”
    - .assay_sc (or .obs_key_assay_sc): Protocol with which data was collected.
        Example: self.assay_sc = “10x”

Optional (if available):

    - .age (or .obs_key_age): Age of individual sampled.
        Example: self.age = 80  # (80 years old for human)
    - .dev_stage (or .obs_key_dev_stage): Developmental stage of individual sampled.
        Example: self.dev_stage = “mature”
    - .ethnicity (or .obs_key_ethnicity): Ethnicity of individual sampled (only for human).
        Example: self.ethnicity = “free text”
    - .healthy (or .obs_key_healthy): Is the sampled from a disease individual? (bool)
        Example: self.healthy = True
    - .sex (or .obs_key_sex): Sex of individual sampled.
        Example: self.sex = “male”
    - .state_exact (or .obs_key_state_exact): Exact disease state
        self.state_exact = free text
    - .obs_key_cellontology_original: Column in .obs in which free text cell type names are stored.
        Example: self.obs_key_cellontology_original = 'CellType'
    - .year: Year of publication:
        Example: self.year = 2019
    - .cell_line: Which cell line was used for the experiment (for cell culture samples)
        Example: self.cell_line = "409B2 (CVCL_K092)"
    - .assay_differentiation: Which protocol was used for the differentiation of the cells (for cell culture samples)
    - .assay_type_differentiation: Which protocol-type was used for the differentiation of the cells: guided or unguided
        (for cell culture samples)

How do I cache data sets?
~~~~~~~~~~~~~~~~~~~~~~~~~
When loading a dataset with `Dataset.load(),`you can specify if the adata object
should be cached or not  (allow_caching= True).
If set to True, the loaded adata object will be cached as an h5ad object for faster reloading.

How do I add cell type annotation?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We are simplifying this right now, new instructions will be available second half of January.

Why are constructor (`__init__`) and loading function (`_load`) split in the template data loader?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Initiation and data set loading are handled separately to allow lazy loading.
All steps that are required to load the count data and
additional metadata should be defined solely in the `_load` section.
Setting of class metadata such as `.doi`, `.id` etc. should be done in the constructor.

How do I tell sfaira where the gene names are?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
By setting the attributes `.var_symbol_col` or `.var_ensembl_col` in the constructor.
If the gene names are in the index of this data frame, you can set “index” as the value of these attributes.

I only have gene symbols (human readable names, often abbreviations), such as HGNC or MGI, but not ENSEMBL identifiers, is that a problem?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No, that is not a problem. They will automatically be converted to Ensembl IDs.
You can, however, specify the reference genome in `Dataset.load(match_to_reference = ReferenceGenomeName)`
to which the names should be mapped to.

I have CITE-seq data, where can I put the protein quantification?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will soon provide a structured interface for loading and accessing CITE-seq data,
for now you can add it into `self.adata.obsm[“CITE”]`.
