Contributing data
==================

Adding datasets to sfaira is a great way to increase the visibility of your dataset and to make it available to a large audience.
This process requires a couple of steps as outlined in the following sections.


.. figure:: https://user-images.githubusercontent.com/21954664/126300611-c5ba18b7-7c88-4bb1-8865-a20587cd5f7b.png
   :alt: sfaira adding datasets

   Overview of contributing dataloaders to sfaira. First, ensure that your data is not yet available as a dataloader.
   Next, create a dataloader and validate it. Afterwards, annotate it to finally test it. Finally, submit your dataloader to sfaira.

sfaira features an interactive way of creating, formatting and testing dataloaders through a command line interface (CLI).
The common workflow using the CLI looks as follows:

1. Check that the data loader was not already implemented.
    We will open issues for all planned data loaders, so you can search both the code_ base and our GitHub issues_ for
    matching data loaders before you start writing one.
    The core data loader identified is the directory compatible doi,
    which is the doi with all special characters replaced by "_" and a "d" prefix is used:
    "10.1016/j.cell.2019.06.029" becomes "d10_1016_j_cell_2019_06_029".
    Searching for this string should yield a match if it is already implemented, take care to look for both
    preprint and publication DOIs if both are available.
    We will also mention publication names in issues, you will however not find these in the code.

.. _code: https://github.com/theislab/sfaira/tree/dev
.. _issues: https://github.com/theislab/sfaira/issues

2. Install sfaira.
    Clone sfaira into a local repository from `dev` branch and install via pip.

.. code-block::

    cd target_directory
    git clone https://github.com/theislab/sfaira.git
    git checkout dev
    # git pull  # use this to update your installation
    cd sfaira  # go into sfaira directory
    pip install -e .  # install
..

3. Create a new dataloader.
    When creating a dataloader with ``sfaira create-dataloader`` dataloader specific attributes such as organ, organism
    and many more are prompted for.
    We provide a description of all meta data items at the bottom of this file.
    If the requested information is not available simply hit enter and continue until done.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    git checkout -b YOUR_BRANCH_NAME  # create a new branch for your data loader.
    sfaira create-dataloader


The created files are created in the sfaira installation under `sfaira/data/dataloaders/loaders/--DOI-folder--`,
where the DOI-specific folder starts with `d` and is followed by the DOI in which all special characters are replaced
by `_`, below referred to as `--DOI-folder--`:

.. code-block::

    ├──sfaira/data/dataloaders/loaders/--DOI-folder--
        ├── extra_description.txt <- Optional extra description file
        ├── __init__.py
        ├── NA_NA_2021_NA_Einstein_001.py <- Contains the load function to load the data
        ├── NA_NA_2021_NA_Einstein_001.yaml <- Specifies all data loader data
..

4. Correct yaml file.
    Correct errors in `sfaira/data/dataloaders/loaders/--DOI-folder--/NA_NA_2021_NA_Einstein_001.yaml` file and add
    further attributes you may have forgotten in step 2.
    This step is optional.

5. Make downloaded data available to sfaira data loader testing.
    Identify the raw files as indicated in the dataloader classes and copy them into your directory structure as
    required by your data loader.
    Note that this should be the exact files that are uploaded to cloud servers such as GEO:
    Do not decompress these files ff these files are archives such as zip, tar or gz.
    Instead, navigate the archives directly in the load function (step 5).
    Copy the data into `sfaira/unit_tests/template_data/--DOI-folder--/`.
    This folder is masked from git and only serves for temporarily using this data for loader testing.
    After finishing loader contribution, you can delete this data again without any consequences for your loader.

6. Write load function.
    Fill load function in `sfaira/data/dataloaders/loaders/--DOI-folder--NA_NA_2021_NA_Einstein_001.py`.

7. Validate the dataloader with the CLI.
    Next validate the integrity of your dataloader content with ``sfaira validate-dataloader <path to *.yaml>``.
    All tests must pass! If any of the tests fail please revisit your dataloader and add the missing information.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    sfaira validate-dataloader <path>``
..

8. Create cell type annotation if your data set is annotated.
    Note that this will abort with error if there are bugs in your data loader.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    # sfaira annotate <path>`` TODO
..

9. Mitigate automated cell type maps.
        Sfaira creates a cell type mapping `.tsv` file in the directory in which your data loaders is located if you
        indicated that annotation is present by filling `cell_types_original_obs_key`.
        This file is: `NA_NA_2021_NA_Einstein_001.tsv`.
        This file contains two columns with one row for each unique cell type label.
        The free text identifiers in the first column "source",
        and the corresponding ontology term in the second column "target".
        You can write this file entirely from scratch.
        Sfaira also allows you to generate a first guess of this file using fuzzy string matching
        which is automatically executed when you run the template data loader unit test for the first time with you new
        loader.
        Conflicts are not resolved in this first guess and you have to manually decide which free text field corresponds
        to which ontology term in the case of conflicts.
        Still, this first guess usually drastically speeds up this annotation harmonization.
        Note that you do not have to include the non-human-readable IDs here as they are added later in a fully
        automated fashion.

10. Test data loader.
        Note that this will abort with error if there are bugs in your data loader.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    # sfaira test-dataloader <path>`` TODO
..

11. Make loader public.
        You can contribute the data loader to public sfaira as code through a pull request.
        Note that you can also just keep the data loader in your local installation or keep it in sfaira_extensions
        if you do not want to make it public.
        Note that we do not manage data upload!
        During publication, you would upload this data set to a server like GEO and the data loader contributed to
        sfaira would use this download link.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    git add *
    git commit  # enter your commit description
    # Next make sure you are up to date with dev
    git checkout dev
    git pull
    git checkout YOUR_BRANCH_NAME
    git merge dev
    git push  # this starts the pull request.
..

The following sections will first describe the underlying design principles of sfaira dataloaders and
then explain how to interactively create, validate and test dataloaders.


Writing dataloaders
---------------------

The study-centric data loader module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the sfaira code, data loaders are organised into directories, which correspond to publications.
All data loaders corresponding to data sets of one study are grouped into this directory.
Next, each data set is represented by one data loader python file in this directory.
See below for more complex set ups with repetitive data loader code.


The data loader python file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each data set, ie a single file or a set of files with similar structures, has its own data loader function and a yaml
files that describes its meta data.
Alternatively to the (preferred) yaml file, meta data can be also be described in a constructor of a class in the same python file
as the loading function. For a documentation on writing a python class-based dataloader, please see here: https://github.com/theislab/sfaira/blob/dev/docs/adding_dataset_classes.rst
A detailed description of all meta data is given at the bottom of this page.

1. A yaml file or constructor of the following form that can be used to interact with the data set
before it is loaded into memory:

.. code-block:: yaml

    dataset_structure:
        dataset_index: 1
        sample_fns:
    dataset_wise:
        author:
        doi:
        download_url_data:
        download_url_meta:
        normalization:
        primary_data:
        year:
    dataset_or_observation_wise:
        assay_sc:
        assay_sc_obs_key:
        assay_differentiation:
        assay_differentiation_obs_key:
        assay_type_differentiation:
        assay_type_differentiation_obs_key:
        bio_sample:
        bio_sample_obs_key:
        cell_line:
        cell_line_obs_key:
        development_stage:
        development_stage_obs_key:
        disease_stage:
        disease_obs_key:
        ethnicity:
        ethnicity_obs_key:
        individual:
        individual_obs_key:
        organ:
        organ_obs_key:
        organism:
        organism_obs_key:
        sample_source:
        sample_source_obs_key:
        sex:
        sex_obs_key:
        state_exact:
        state_exact_obs_key:
        tech_sample:
        tech_sample_obs_key:
    observation_wise:
        cell_types_original_obs_key:
    feature_wise:
        gene_id_ensembl_var_key:
        gene_id_symbols_var_key:
    meta:
        version: "1.0"


2. A function called to load the data set into memory:
It is important to set an automated path indicating the location of the raw files here.
Our recommendation for this directory set-up is that you define a directory folder in your directory structure
in which all of these raw files will be (self.path) and then add a sub-directory named as
`self.directory_formatted_doi` (ie. the doi with all special characters replaced by "_" and place the raw files
directly into this sub directory.

.. code-block:: python

    def load(data_dir, fn=None) -> anndata.AnnData:
        fn = os.path.join(data_dir, "my.h5ad")
        adata = anndata.read(fn)  # loading instruction into adata, use other ones if the data is not h5ad
        return adata

In summary, a the dataloader for a mouse lung data set could look like this:

.. code-block:: yaml

    dataset_structure:
        dataset_index: 1
        sample_fns:
    dataset_wise:
        author: "me"
        doi:
            - "my preprint"
            - "my peer-reviewed publication"
        download_url_data: "my GEO upload"
        download_url_meta:
        normalization: "raw"
        primary_data:
        year:
    dataset_or_observation_wise:
        assay_sc: "smart-seq2"
        assay_sc_obs_key:
        assay_differentiation:
        assay_differentiation_obs_key:
        assay_type_differentiation:
        assay_type_differentiation_obs_key:
        bio_sample:
        bio_sample_obs_key:
        cell_line:
        cell_line_obs_key:
        development_stage:
        development_stage_obs_key:
        disease_stage:
        disease_obs_key:
        ethnicity:
        ethnicity_obs_key:
        individual:
        individual_obs_key:
        organ: "lung"
        organ_obs_key:
        organism: "mouse"
        organism_obs_key:
        sample_source: "primary_tissue"
        sample_source_obs_key:
        sex:
        sex_obs_key:
        state_exact:
        state_exact_obs_key:
        tech_sample:
        tech_sample_obs_key:
    observation_wise:
        cell_types_original_obs_key: "louvain_named"
    feature_wise:
        gene_id_ensembl_var_key:
        gene_id_symbols_var_key:
    meta:
        version: "1.0"

.. code-block:: python

    def load(data_dir, fn=None) -> anndata.AnnData:
        fn = os.path.join(data_dir, "my.h5ad")
        adata = anndata.read(fn)
        return adata


Data loaders can be added into a copy of the sfaira repository and can be used locally before they are contributed to
the public sfaira repository.
Alternatively, we also provide the optional dependency sfaira_extensions (https://github.com/theislab/sfaira_extension)
in which local data and cell type annotation can be managed separately but still be loaded as usual through sfaira.
The data loaders and cell type annotation formats between sfaira and sfaira_extensions are identical and can be easily
copied over.

Loading multiple files of similar structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only one loader has to be written for each set of files that are similarly structured which belong to one DOI.
`sample_fns` in `dataset_structure` in the `.yaml` indicates the presence of these files.
The identifiers listed there do not have to be the full file names.
They are received by `load()`  as the argument `sample_fn` and can then be used in custom code in `load()` to load
the correct file.
This allows sharing code across these files in `load()`.
If these files share all meta data in the `.yaml`, you do not have to change anything else here.
If a some meta data items are file specific, you can further subdefine them under the keys in this `.yaml` via their
identifiers stated here.
In the following example, we show how this formalism can be used to identify one file declared as "A" as a healthy
lung sample and another file "B" as a healthy pancreas sample.

.. code-block:: python

    dataset_structure:
        dataset_index: 1
        sample_fns:
            - "A"
            - "B"
    dataset_wise:
        # ... part of yaml omitted ...
    dataset_or_observation_wise:
        # ... part of yaml omitted
        healthy: True
        healthy_obs_key:
        individual:
        individual_obs_key:
        organ:
            A: "lung"
            B: "pancreas"
        organ_obs_key:
        # part of yaml omitted ...
..

Note that not all meta data items have to subdefined into "A" and "B" but only the ones with differing values!
The corresponding `load` function would be:

.. code-block:: python

    def load(data_dir, sample_fn, fn=None) -> anndata.AnnData:
        # The following reads either my_file_A.h5ad or my_file_B.h5ad which correspond to A and B in the yaml.
        fn = os.path.join(data_dir, f"my_file_{sample_fn}.h5ad")
        adata = anndata.read(fn)
        return adata
..


Loading third party annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, the data set in question is already in the sfaira zoo but there is alternative (third party), cell-wise
annotation of the data.
This could be different cell type annotation for example.
The underlying data (count matrix and variable names) stay the same in these cases, and often, even some cell-wise
meta data are kept and only some are added or replaced.
Therefore, these cases do not require an additional `load()` function.
Instead, you can contribute `load_annotation_*()` functions into the `.py` file of the corresponding study.
You can chose an arbitrary suffix for the function but ideally one that identifies the source of this additional
annotation in a human readable manner at least to someone who is familiar with this data set.
Second you need to add this function into the dictionary `LOAD_ANNOTATION` in the `.py` file, with the suffix as a key.
If this dictionary does not exist yet, you need to add it into the `.py` file with this function as its sole entry.
Here an example of a `.py` file with additional annotation:

.. code-block:: python

    def load(data_dir, sample_fn, **kwargs):
        pass

    def load_annotation_meta_study_x(data_dir, sample_fn, **kwargs):
        # Read a tabular file indexed with the observation names used in the adata used in load().
        pass

    def load_annotation_meta_study_y(data_dir, sample_fn, **kwargs):
        # Read a tabular file indexed with the observation names used in the adata used in load().
        pass

    LOAD_ANNOTATION = {
        "meta_study_x": load_annotation_meta_study_x,
        "meta_study_y": load_annotation_meta_study_y,
    }


The table returned by `load_annotation_meta_study_x` needs to be indexed with the observation names used in `.adata`,
the object generated in `load()`.
If `load_annotation_meta_study_x` contains a subset of the observations defined in `load()`,
and this alternative annotation is chosen,
`.adata` is subsetted to these observations during loading.

You can also add functions in the `.py` file in the same DOI-based module in sfaira_extensions if you want to keep this
additional annotation private.
For this to work with a public data loader, you need nothing more than the `.py` file with this `load_annotation_*()`
function and the `LOAD_ANNOTATION` of these private functions in sfaira_extensions.

To access additional annotation during loading, use the setter functions `additional_annotation_key` on an instance of
either `Dataset`, `DatasetGroup` or `DatasetSuperGroup` to define data sets
for which you want to load additional annotation and which additional you want to load for these.
See also the docstrings of these functions for further details on how these can be set.


Creating dataloaders with the commandline interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sfaira features an interactive way of creating, formatting and testing dataloaders.
The common workflow look as follows:

1. Create a new dataloader with ``sfaira create-dataloader``
2. Validate the dataloader with ``sfaira lint-dataloader <path>``
3. Test the dataloader with ``sfaira test-dataloader . --doi <doi> --test-data <folder_above_test_data>``

When creating a dataloader with ``sfaira create-dataloader`` common information such as
your name and email are prompted for, followed by dataloader specific attributes such as organ, organism and many more.
If the requested information is not available simply hit enter and continue until done. If you have mixed organ or organism
data you will have to resolve this manually later. Your dataloader template will be created in your current working directory
in a folder resembling your doi.

The created files are:

.. code-block::

    ├── extra_description.txt <- Optional extra description file
    ├── __init__.py
    ├── NA_NA_2021_NA_Einstein_001.py <- Contains the load function to load the data
    ├── NA_NA_2021_NA_Einstein_001.yaml <- Specifies all data loader data

Now simply fill in all missing properties in your dataloader scripts and yaml file.
When done optionally run ``sfaira clean-dataloader <path to *.yaml>`` on the just filled out dataloader yaml file.
All unused attributes will be removed.

Next validate the integrity of your dataloader content with ``sfaira validate-dataloader <path to *.yaml>``.
All tests must pass! If any of the tests fail please revisit your dataloader and add the missing information.

Finally, copy your dataloader into the ``sfaira/dataloaders/loaders/`` folder.
Now you can test your dataloader with ``sfaira test-dataloader <path_to_sfaira> --doi <doi> --test-data <template_data_folder>``.
Note that sfaira expects a folder structure for the test data such as:

.. code-block::

    ├── template_data
    │   └── d10_1016_j_cmet_2019_01_021
    │       ├── GSE117770_RAW.tar
    │       ├── GSM3308545_NOD_08w_A_annotation.csv
    │       ├── GSM3308547_NOD_08w_C_annotation.csv
    │       ├── GSM3308548_NOD_14w_A_annotation.csv
    │       ├── GSM3308549_NOD_14w_B_annotation.csv
    │       ├── GSM3308550_NOD_14w_C_annotation.csv
    │       ├── GSM3308551_NOD_16w_A_annotation.csv
    │       ├── GSM3308552_NOD_16w_B_annotation.csv
    │       └── GSM3308553_NOD_16w_C_annotation.csv

Pass the path to the template_data folder, not the doi. Sfaira will use this path to cache further data for speedups.
All tests must pass! If any of the tests fail please revisit your dataloader and fix the error.

Map cell type labels to ontology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The entries in `self.cell_types_original_obs_key` are free text but are mapped to an ontology via a .tsv file with
the same name and directory as the python file in which the data loader is located.
This .tsv contains two columns with one row for each unique cell type label.
The free text identifiers in the first column "source",
and the corresponding ontology term in the second column "target".
You can write this file entirely from scratch.
Sfaira also allows you to generate a first guess of this file using fuzzy string matching
which is automatically executed when you run the template data loader unit test for the first time with you new loader.
Conflicts are not resolved in this first guess and you have to manually decide which free text field corresponds to which
ontology term in the case of conflicts.
Still, this first guess usually drastically speeds up this annotation harmonization.

Cell type ontology management
-----------------------------

Sfaira maintains a wrapper of the Cell Ontology as a class which allows additions to this ontology.
This allows us to use the core ontology used in the community as a backbone and to keep up with newly identifed cell types on our own.
We require all extensions of the core ontology not to break the directed acyclic graph that is the ontology:
Usually, such extensions would be additional leave nodes.

Second, we maintain cell type universes for anatomic structures.
These are dedicated for cell type-dependent models which require a defined set of cell types.
Such a universe is a set of nodes in the ontology.

Contribute cell types to ontology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please open an issue on the sfaira repo with a description what type of cell type you want to add.


Metadata
--------

Required fields
~~~~~~~~~~~~~~~

Most meta data fields are optional in sfaira.
Required are:

- dataset_structure: dataset_index is required.
- dataset_wise: author, doi, download_url_data, normalisation and year are required.
- dataset_or_observation_wise: organism is required.
- observation_wise: None are required.
- feature_wise: gene_id_ensembl_var_key or gene_id_symbols_var_key is required.
- misc: None are required.

Field descriptions
~~~~~~~~~~~~~~~~~~

We constrain meta data by ontologies where possible.
Meta data can either be dataset-wise, observation-wise or feature-wise.

Dataset structure meta data are in the section `dataset_structure` in the `.yaml` file.

- dataset_index [int]
    Numeric identifier of the first loader defined by this python file.
    Only relevant if multiple python files for one DOI generate loaders of the same name.
    In these cases, this numeric index can be used to distinguish them.
- sample_fns [list of strings]
    If there are multiple data files which can be covered by one `load()` function and `.yaml` file because they are
    structured similarly, these can identified here.
    See also section `Loading multiple files of similar structure`.

Dataset-wise meta data are in the section `dataset_wise` in the `.yaml` file.

- author [list of strings]
    List of author names of dataset (not of loader).
- doi [list of strings]
    DOIs associated with dataset.
    These can be preprints and journal publication DOIs.
- download_url_data [list of strings]
    Download links for data.
    Full URLs of all data files such as count matrices. Note that distinct observation-wise annotation files can be
    supplied in download_url_meta.
- download_url_meta [list of strings]
    Download links for observation-wise data.
    Full URLs of all observation-wise meta data files such as count matrices.
    This attribute is optional and not necessary ff observation-wise meta data is already in the files defined in
    `download_url_data`, e.g. often the case for .h5ad`.
- normalization: Data normalisation {"raw", "scaled"}
    Type of normalisation of data stored in `adata.X` emitted by the `load()` function.
- year: Year in which sample was first described [integer]
    Pre-print publication year.

Meta-data which can either be dataset- or observation-wise are in the section `dataset_or_observation_wise` in the
`.yaml` file.
They can all be supplied as `NAME` or as `NAME_obs_key`:
The former indicates that the entire data set has the value stated in the yaml.
The latter, `NAME_obs_key`, indicates that there is a column in `adata.obs` emitted by the `load()` function of the name
`NAME_obs_key` which contains the annotation per observation for this meta data item.
Note that in both cases the value, or the column values, have to fulfill contraints imposed on the meta data item as
outlined below.

- assay_sc and assay_sc_obs_key [ontology term]
    Choose a term from https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0010183&viewMode=All&siblings=false
- assay_differentiation and assay_differentiation_obs_key [string]
    Try to provide a base differentiation protocol (eg. "Lancaster, 2014") as well as any amendments to the original
    protocol.
- assay_type_differentiation and assay_type_differentiation_obs_key {"guided", "unguided"}
    For cell-culture samples: Whether a guided (patterned) differentiation protocol was used in the experiment.
- bio_sample and bio_sample_obs_key [string]
    Column name in `adata.obs` emitted by the `load()` function which reflects biologically distinct samples, either
    different in condition or biological replicates, as a categorical variable.
    The values of this column are not constrained and can be arbitrary identifiers of observation groups.
    You can concatenate multiple columns to build more fine grained observation groupings by concatenating the column
    keys with `*` in this string, e.g. `patient*treatment` to get one `bio_sample` for each patient and treatment.
    Note that the notion of biologically distinct sample is slightly subjective, we allow this element to allow
    researchers to distinguish technical and biological replicates within one study for example.
    See also the meta data items `individual` and `tech_sample`.
- cell_line and cell_line_obs_key [ontology term]
    Cell line name from the cellosaurus cell line database (https://web.expasy.org/cellosaurus/)
- developmental_stage and developmental_stage_obs_key [ontology term]
    Developmental stage (age) of individual sampled.
    Choose from HSAPDV (https://www.ebi.ac.uk/ols/ontologies/hsapdv) for human
    or from MMUSDEV (https://www.ebi.ac.uk/ols/ontologies/mmusdv) for mouse.
- disease and disease_obs_key [ontology term]
    Choose from MONDO (https://www.ebi.ac.uk/ols/ontologies/mondo) for human
- ethnicity and ethnicity_obs_key [ontology term]
    Choose from HANCESTRO (https://www.ebi.ac.uk/ols/ontologies/hancestro)
- individual and individual_obs_key [string]
    Column name in `adata.obs` emitted by the `load()` function which reflects the indvidual sampled as a categorical
    variable.
    The values of this column are not constrained and can be arbitrary identifiers of observation groups.
    You can concatenate multiple columns to build more fine grained observation groupings by concatenating the column
    keys with `*` in this string, e.g. `group1*group2` to get one `individual` for each group1 and group2 entry.
    Note that the notion of individuals is slightly mal-defined in some cases, we allow this element to allow
    researchers to distinguish sample groups that originate from biological material with distinct genotypes.
    See also the meta data items `individual` and `tech_sample`.
- organ and organ_obs_key [ontology term]
    The UBERON anatomic location of the sample (https://www.ebi.ac.uk/ols/ontologies/uberon).
- organism and organism_obs_key. {"mouse", "human"}.
    The organism from which the sample originates.
    In the future, we will use NCBITAXON (https://www.ebi.ac.uk/ols/ontologies/ncbitaxon).
- primary_data [bool]
    Whether contains cells that were measured in this study (ie this is not a meta study on published data).
- sample_source and sample_source_obs_key. {"primary_tissue", "2d_culture", "3d_culture", "tumor"}
    Which cellular system the sample was derived from.
- sex and sex_obs_key. Sex of individual sampled. {"female", "male", None}
    Sex of the individual sampled.
- state_exact and state_exact_obs_key [string]
    Free text description of condition.
    If you give treatment concentrations, intervals or similar measurements use square brackets around the quantity
    and use units: `[1g]`
- tech_sample and tech_sample_obs_key [string]
    Column name in `adata.obs` emitted by the `load()` function which reflects technically distinct samples, either
    different in condition or technical replicates, as a categorical variable.
    Any data batch is a `tech_sample`.
    The values of this column are not constrained and can be arbitrary identifiers of observation groups.
    You can concatenate multiple columns to build more fine grained observation groupings by concatenating the column
    keys with `*` in this string, e.g. `patient*treatment*protocol` to get one `tech_sample` for each patient, treatment
    and measurement protocol.
    See also the meta data items `individual` and `tech_sample`.

Meta-data which are strictly observation-wise are in the section `observation_wise` in the `.yaml` file:

- cell_types_original_obs_key [string]
    Column name in `adata.obs` emitted by the `load()` function which contains free text cell type labels.

Meta-data which are feature-wise are in the section `feature_wise` in the `.yaml` file:

- gene_id_ensembl_var_key [string]
    Name of the column in `adata.var` emitted by the `load()` which contains ENSEMBL gene IDs.
    This can also be "index" if the ENSEMBL gene names are in the index of the `adata.var` data frame.
- gene_id_symbols_var_key:.[string]
    Name of the column in `adata.var` emitted by the `load()` which contains gene symbol:
    HGNC for human and MGI for mouse.
    This can also be "index" if the gene symbol are in the index of the `adata.var` data frame.

The meta data on the meta data file do not have to modified by you are automatically controlled are in the section
`meta` in the `.yaml` file:

- version: [string]
    Version identifier of meta data scheme.

The class-based data loader python file
----------------------------------------
As an alternative to the preferred yaml-based dataloaders, users can provide a dataloader class together with the load function.
In this scenario, meta data is described in a constructor of a class in the same python file as the loading function.

1. A constructor of the following form that contains all the relevant metadata that is available before the actual dataset is loaded to memory.

.. code-block:: python

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        # Data set meta data: You do not have to include all of these and can simply skip lines corresponding
        # to attritbutes that you do not have access to. These are meta data on a sample level.
        # The meta data attributes labeled with (*) may als be supplied per cell, see below,
        # in this case, if you supply a .obs_key* attribute, you ccan leave out the sample-wise attribute.

        self.id = x  # unique identifier of data set (Organism_Organ_Year_AssaySc_NumberOfDataset_FirstAuthorLastname_doi).

        self.author = x  # author (list) who sampled / created the data set
        self.doi = x  # doi of data set accompanying manuscript

        self.download_url_data = x  # download website(s) of data files
        self.download_url_meta = x  # download website(s) of meta data files

        self.assay_sc = x  # (*, optional) protocol used to sample data (e.g. smart-seq2)
        self.assay_differentiation = x  # (*, optional) protocol used to differentiate the cell line (e.g. Lancaster, 2014)
        self.assay_type_differentiation = x  # (*, optional) type of protocol used to differentiate the cell line (guided/unguided)
        self.cell_line = x # (*, optional) cell line used (for cell culture samples)
        self.dev_stage = x  # (*, optional) developmental stage of organism
        self.ethnicity = x  # (*, optional) ethnicity of sample
        self.healthy = x  # (*, optional) whether sample represents a healthy organism
        self.normalisation = x  # (optional) normalisation applied to raw data loaded (ideally counts, "raw")
        self.organ = x  # (*, optional) organ (anatomical structure)
        self.organism = x  # (*) species / organism
        self.sample_source = x  # (*) whether the sample came from primary tissue or cell culture
        self.sex = x  # (*, optional) sex
        self.state_exact = x  # (*, optional) exact disease, treatment or perturbation state of sample
        self.year = x  # year in which sample was acquired

        # The following meta data may instead also be supplied on a cell level if an appropriate column is present in the
        # anndata instance (specifically in .obs) after loading.
        # You need to make sure this is loaded in the loading script)!
        # See above for a description what these meta data attributes mean.
        # Again, if these attributes are note available, you can simply leave this out.
        self.obs_key_assay_sc = x  # (optional, see above, do not provide if .assay_sc is provided)
        self.obs_key_assay_differentiation = x  # (optional, see above, do not provide if .age is assay_differentiation)
        self.obs_key_assay_type_differentiation = x  # (optional, see above, do not provide if .assay_type_differentiation is provided)
        self.obs_key_cell_line = x # (optional, see above, do not provide if .cell_line is provided)
        self.obs_key_dev_stage = x  # (optional, see above, do not provide if .dev_stage is provided)
        self.obs_key_ethnicity = x  # (optional, see above, do not provide if .ethnicity is provided)
        self.obs_key_healthy = x  # (optional, see above, do not provide if .healthy is provided)
        self.obs_key_organ = x  # (optional, see above, do not provide if .organ is provided)
        self.obs_key_organism = x  # (optional, see above, do not provide if .organism is provided)
        self.obs_key_sample_source = x  # (optional, see above, do not provide if .sample_source is provided)
        self.obs_key_sex = x  # (optional, see above, do not provide if .sex is provided)
        self.obs_key_state_exact = x  # (optional, see above, do not provide if .state_exact is provided)
        # Additionally, cell type annotation is ALWAYS provided per cell in .obs, this annotation is optional though.
        # name of column which contain streamlined cell ontology cell type classes:
        self.obs_key_cell_types_original = x  # (optional)
        # This cell type annotation is free text but is mapped to an ontology via a .tsv file with the same name and
        # directory as the python file of this data loader (see below).


2. A function called to load the data set into memory:
It is important to set an automated path indicating the location of the raw files here.
Our recommendation for this directory set-up is that you define a directory folder in your directory structure
in which all of these raw files will be (self.path) and then add a sub-directory named as
`self.directory_formatted_doi` (ie. the doi with all special characters replaced by "_" and place the raw files
directly into this sub directory.

.. code-block:: python

    def load(data_dir, fn=None) -> anndata.AnnData:
        fn = os.path.join(data_dir, "my.h5ad")
        adata = anndata.read(fn)  # loading instruction into adata, use other ones if the data is not h5ad
        return adata

In summary, a python file for a mouse lung data set could look like this:

.. code-block:: python

    class MyDataset(DatasetBase)
        def __init__(
                self,
                path: Union[str, None] = None,
                meta_path: Union[str, None] = None,
                cache_path: Union[str, None] = None,
                **kwargs
        ):
            super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
            self.author = "me"
            self.doi = ["my preprint", "my peer-reviewed publication"]
            self.download_url_data = "my GEO upload"
            self.normalisation = "raw"  # because I uploaded raw counts, which is good practice!
            self.organ = "lung"
            self.organism = "mouse"
            self.assay_sc = "smart-seq2"
            self.year = "2020"
            self.sample_source = "primary_tissue"

            self.obs_key_cell_types_original = "louvain_named"  # i save my cell type names in here

    def load(data_dir, fn=None) -> anndata.AnnData:
        fn = os.path.join(data_dir, "my.h5ad")
        adata = anndata.read(fn)
        return adata
