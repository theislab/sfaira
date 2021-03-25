Adding data sets
===================

Adding datasets to sfaira is a great way to increase the visibility of your dataset and to make it available to a large audience.
This process requires a couple of steps as outlined in the following sections.
sfaira features an interactive way of creating, formatting and testing dataloaders through a command line interface (CLI).
The common workflow using the CLI looks as follows:

1. Install sfaira.
    Clone sfaira into a local repository from `dev` branch and install via pip.

.. code-block::

    cd target_directory
    git clone https://github.com/theislab/sfaira.git
    git checkout dev
    # git pull  # use this to update your installation
    cd sfaira  # go into sfaira directory
    pip install -e .  # install
..

2. Create a new dataloader.
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

3. Correct yaml file.
    Correct errors in `sfaira/data/dataloaders/loaders/--DOI-folder--/NA_NA_2021_NA_Einstein_001.yaml` file and add
    further attributes you may have forgotten in step 2.
    This step is optional.

4. Make downloaded data available to sfaira data loader testing.
    Identify the raw files as indicated in the dataloader classes and copy them into your directory structure as
    required by your data loader.
    Note that this should be the exact files that are uploaded to cloud servers such as GEO:
    Do not decompress these files ff these files are archives such as zip, tar or gz.
    Instead, navigate the archives directly in the load function (step 5).
    Copy the data into `sfaira/unit_tests/template_data/--DOI-folder--/`.
    This folder is masked from git and only serves for temporarily using this data for loader testing.
    After finishing loader contribution, you can delete this data again without any consequences for your loader.

5. Write load function.
    Fill load function in `sfaira/data/dataloaders/loaders/--DOI-folder--NA_NA_2021_NA_Einstein_001.py`.

6. Clean the dataloader with a supervicial check (lint).
    This step is optional.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    sfaira clean-dataloader <path to *.yaml>
..

7. Validate the dataloader with the CLI.
    Next validate the integrity of your dataloader content with ``sfaira lint-dataloader <path to *.yaml>``.
    All tests must pass! If any of the tests fail please revisit your dataloader and add the missing information.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    sfaira lint-dataloader <path>``
..

8. Create cell type annotation if your data set is annotated.
    Note that this will abort with error if there are bugs in your data loader.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    # sfaira annotate <path>`` TODO
..

9. Mitigate automated cell type maps.
    Sfaira creates a cell type mapping `.tsv` file in the directory in which your data loaders is located if you
    indicated that annotation is present by filling `cellontology_original_obs_key`.
    This file is: `NA_NA_2021_NA_Einstein_001.tsv`.
    This file contains two columns with one row for each unique cell type label.
    The free text identifiers in the first column "source",
    and the corresponding ontology term in the second column "target".
    You can write this file entirely from scratch.
    Sfaira also allows you to generate a first guess of this file using fuzzy string matching
    which is automatically executed when you run the template data loader unit test for the first time with you new loader.
    Conflicts are not resolved in this first guess and you have to manually decide which free text field corresponds to which
    ontology term in the case of conflicts.
    Still, this first guess usually drastically speeds up this annotation harmonization.
    Note that you do not have to include the non-human-readable IDs here as they are added later in a fully automated
    fashion.

10. Test data loader
    Note that this will abort with error if there are bugs in your data loader.

.. code-block::

    # make sure you are in the top-level sfaira directory from step 1
    # sfaira test <path>`` TODO
..

11. You can contribute the data loader to public sfaira as code through a pull request.
    Note that you can also just keep the data loader in your local installation or keep it in sfaira_extensions if you
    do not want to make it public.
    Note that we do not manage data upload!
    During publication, you would upload this data set to a server like GEO and the data loader contributed to sfaira
    would use this download link.

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

Check that the data loader was not already implemented
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will open issues for all planned data loaders, so you can search both the code_ base and our GitHub issues_ for
matching data loaders before you start writing one.
The core data loader identified is the directory compatible doi,
which is the doi with all special characters replaced by "_" and a "d" prefix is used:
"10.1016/j.cell.2019.06.029" becomes "d10_1016_j_cell_2019_06_029".
Searching for this string should yield a match if it is already implemented, take care to look for both
preprint and publication DOIs if both are available. We will also mention publication names in issues, you will however not find these in the code.

.. _code: https://github.com/theislab/sfaira/tree/dev
.. _issues: https://github.com/theislab/sfaira/issues


The data loader python file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each data set, ie a single file or a set of files with similar structures, has its own data loader function and a yaml
files that describes its meta data.
Alternatively to the (preffered) yaml file, meta data can be also be described in a constructor of a class in the same python file
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
        year:
    dataset_or_observation_wise:
        age:
        age_obs_key:
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
        healthy:
        healthy_obs_key:
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
        cellontology_original_obs_key:
    feature_wise:
        var_ensembl_col:
        var_symbol_col:
    misc:
        healthy_state_healthy:
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
        year:
    dataset_or_observation_wise:
        age:
        age_obs_key:
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
        healthy:
        healthy_obs_key:
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
        cellontology_original_obs_key: "louvain_named"
    feature_wise:
        var_ensembl_col:
        var_symbol_col:
    misc:
        healthy_state_healthy:
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

Using ontologies to train cell type classifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cell type classifiers can be trained on data sets with different coarsity of cell type annotation using aggregate
cross-entropy as a loss and aggregate accuracy as a metric.
The one-hot encoded cell type label matrix is accordingly modified in the estimator class in data loading if terms
that correspond to intermediate nodes (rather than leave nodes) are encountered in the label set.

Metadata management
-------------------

We constrain meta data by ontologies where possible. The current restrictions are:

    - .age: unconstrained string
        Use
            - units of years for humans,
            - the E{day} nomenclature for mouse embryos
            - the P{day} nomenclature for young post-natal mice
            - units of weeks for mice older than one week and
            - units of days for cell culture samples.
    - .assay_sc: EFO-constrained string
        Choose a term from https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0010183&viewMode=All&siblings=false
    - .assay_differentiation: unconstrained string
        Try to provide a base differentiation protocol (eg. "Lancaster, 2014") as well as any amendments to the original protocol.
    - .assay_type_differentiation: constrained string, {"guided", "unguided"}
        For cell-culture samples: Whether a guided (patterned) differentiation protocol was used in the experiment.
    - .developmental_stage: unconstrained string
        This will constrained to an ontology in the future,
        try choosing from HSAPDV (https://www.ebi.ac.uk/ols/ontologies/hsapdv) for human
        or from MMUSDEV (https://www.ebi.ac.uk/ols/ontologies/mmusdv) for mouse.
    - .cell_line: cellosaurus-constrained string
        Cell line name from the cellosaurus cell line database (https://web.expasy.org/cellosaurus/)
    - .ethnicity: unconstrained string, this will constrained to an ontology in the future.
        Try choosing from HANCESTRO (https://www.ebi.ac.uk/ols/ontologies/hancestro)
    - .healthy: bool
        Whether the sample is from healthy tissue ({True, False}).
    - .normalisation: unconstrained string, this will constrained to an ontology in the future,
        Try to use {"raw", "scaled"}.
    - .organ: UBERON-constrained string
        The anatomic location of the sample (https://www.ebi.ac.uk/ols/ontologies/uberon).
    - .organism: constrained string, {"mouse", "human"}.
        The organism from which the sample originates.
        In the future, we will use NCBITAXON (https://www.ebi.ac.uk/ols/ontologies/ncbitaxon).
    - .sample_source: constrained string, {"primary_tissue", "2d_culture", "3d_culture", "tumor"}
        Which cellular system the sample was derived from.
    - .sex: constrained string, {"female", "male", None}
        Sex of the individual sampled.
    - .state_exact: unconstrained string, try to be concise and anticipate that this field is queried by automatised searches.
        If you give treatment concentrations, intervals or similar measurements use square brackets around the quantity
        and use units: `[1g]`
    - .year: must be an integer year, e.g. 2020
        Year in which sample was first described (e.g. pre-print publication).

Follow this issue_ for details on upcoming ontology integrations.

.. _issue: https://github.com/theislab/sfaira/issues/16
