.. _adding_data_rst:

Writing data loaders
=====================

For a high-level overview of data management in sfaira, read :ref:`data_life_cycle_rst` first.
In brief, a data loader is a set of instructions allow for streamlining of raw count matrices
and meta data into objects of a target format.
Here, streamlining means that gene names are controlled based on a genome assembly,
metadata items are contrained to follow ontologies,
and key study metadata are described.
This streamlining increases accessibility and visibility of a dataset and to makes it available to a large audience.
In sfaira, data loaders are grouped by scientific study (DOI of a preprint or DOI of a publication).
A data loader for a study is a directory named after the DOI of the study that contains code and text files.
This directory is part of the sfaira python package and, thus, maintained on GitHub.
This allows for data loaders to be maintained via GitHub workflows: contribution and fixes via pull requests and
deployment via repository cloning and package installation.

A dataloader consists of four file components within a single directory:

1. `__init__.py` file which is has same content in all loaders,
2. `ID.py` file that contains a `load()` functions with based instructions of loading raw data on disk,
3. `ID.yaml` file that describes most meta data,
4. `ID*.tsv` files with ontology-wise maps of free-text metadata items to contrained vocabulary.

Note that all dataset-specific components receive and `ID` that is set during the curation process.
Below, we desribe how multiple datasets within a study can be handled with the same dataloder.
In cases where this is not efficient, one can go through the data loader creation process once for each dataset
and then group the resuling loaders (file groups 1-4) in a single directory named after the study's DOI.

An experienced curator can directly write such a data loader.
However, first-time contributors often struggle with the interplay of individual files,
metadata maps from free-text annotation are notoriously buggy
and comprehensive testing is important also for contributions by experienced curators.
Therefore, we broke the process of writing a loader down into phases
and built a CLI to guide users through this process.
Each phase corresponds to one command (one execution of a shell command) in the CLI.
In addition, the CLI guides the user through manual steps that are necessary in each phase.
We structured the process of curation into four phases,
a preparatory phase P precedes CLI execution and is described in this documentation.

- Phase P (``prepare``): data and python environment setup for curation.
- Phase 1 (``create``): a `load()` function (in a `.py`) and a YAML are written.
- Phase 2 (``annotate``): ontology-specific maps of free-text metadata to contrained vocabulary (in `*.tsv`) are written.
- Phase 3 (``finalize``): the data loader is tested and metadata are cleaned up.
- Phase 4 (``upload``): the data loader is uploaded to the sfaira GitHub repository.

An experienced curator could skip using the CLI for phase 1 and write the `.py` and `.yaml` by hand.
In this case, we still highly recommend using the CLI for phase 2 and 3.
Note that phase 2 is only necessary if you have free-text metadata that needs to be mapped,
the CLI will point this out accordingly in phase 1.
This 4-phase cycle completes initial curation and results in data loader code that can be pushed to the sfaira
GitHub repository.
This cycle can be complemented by an optional workflow to cache curated `.h5ad` objects (e.g. on the cellxgene website):

- Phase 5 (``export``): the data loader is used to create a streamlined `.h5ad` of a particular format.
- Phase 6 (``validate-h5ad``): the `.h5ad` from phase 4 is checked for compliance with a particular (e.g. the cellxgene format).

The resuling `.h5ad` can be shared with collaborators or uploaded to data submission servers.

Create a new data loader
-------------------------

Phase P: Preparation
~~~~~~~~~~~~~~~~~~~~~

Before you start writing the data loader, we recommend completing this checks and preparation measures.
Phase P is sub-structured into 3 sub-phases:

* Pa: Name the data loader.
* Pb: Check that the data loader was not already implemented.
* Pc: Prepare an installation of sfaira to use for data loader writing.
* Pd: Download the raw data into a local directory.

Pa. Name the data loader.
    We will decide for a  name of the dataloader based on its DOI.
    Prefix the DOI with `"d"` and replace the special characters in the DOI with `"_"` here to prevent copy mistakes,
    e.g. the DOI `10.1000/j.journal.2021.01.001` becomes `d10_1000_j_journal_2021_01_001`
    Remember to replace this DOI with the DOI of the study you want to contribute, choose a publication (journal)
    DOI if available, otherwise a preprint DOI.
    If neither DOI is available, because this is unpublished data, for example, use an identifier that makes sense to
    you, that is prefixed with `dno_doi` and contains a name of an author of the dataset, e.g.
    `dno_doi_einstein_brain_atlas`.
    We will refer to this name as `DOI-name` and it will be used to label the contributed code and the stored data.

Pb. Check that the data loader was not already implemented.
    We will open issues for all planned data loaders, so you can search both the code_ base and our GitHub issues_ for
    matching data loaders before you start writing one.
    The core data loader identified is the directory compatible doi,
    which is the doi with all special characters replaced by "_" and a "d" prefix is used:
    "10.1016/j.cell.2019.06.029" becomes "d10_1016_j_cell_2019_06_029".
    Searching for this string should yield a match if it is already implemented, take care to look for both
    preprint and publication DOIs if both are available.
    We will also mention publication names in issues, you will however not find these in the code.

Pc. Prepare an installation of sfaira to use for data loader writing.
    Jump to 2d) if you do not require explanations of specifc parts of the shell script.

    1. Install sfaira.
        Clone sfaira into a local repository `DIR_SFAIRA`.

        .. code-block::

            cd DIR_SFAIRA
            git clone https://github.com/theislab/sfaira.git
            cd sfaira
            git checkout dev
        ..
    2. Prepare a local branch of sfaira dedicated to your loader.
        You can name this branch after the `DOI-name`, prefix this branch with `data/` as the code change suggested
        is a data addition.

        .. code-block::

            cd DIR_SFAIRA
            cd sfaira
            git checkout dev
            git pull
            git checkout -b data/DOI-name
        ..
    3. Install sfaira into a conda environment.
        You can for example use pip inside of a conda environment dedicated to data curation.

        .. code-block::

            cd DIR_SFAIRA
            cd sfaira
            git checkout -b data/DOI-name
            conda create -n sfaira_loader
            conda install -n sfaira_loader python=3.8
            conda activate sfaira_loader
            pip install -e .
        ..
    4. Summary of step 1-3.
        P2a-c are all covered by the following code block, remember to name the git branch after your DOI:

        .. code-block::

            cd DIR_SFAIRA
            git clone https://github.com/theislab/sfaira.git
            cd sfaira
            git checkout dev
            git pull
            git checkout -b data/DOI-name
            conda create -n sfaira_loader
            conda install -n sfaira_loader python=3.8
            conda activate sfaira_loader
            pip install -e .
        ..

Pd. Download the raw data into a local directory.
    You will need to set a path in which the data files can be accessed by sfaira, in the following referred to as
    `<path_data>/<DOI-name>/`.
    Identify the raw data files and copy them into the datafolder `<path_data>/<DOI-name>/`.
    Note that this should be the exact files that are downloadable from the download URL you provided in the dataloader:
    Do not decompress these files if these files are archives such as zip, tar or gz.

.. _code: https://github.com/theislab/sfaira/tree/dev/sfaira/data/dataloaders/loaders
.. _issues: https://github.com/theislab/sfaira/issues

Phase 1: create
~~~~~~~~~~~~~~~~

Phase 1 is sub-structured into 2 sub-phases:

* 1a: Create template files (``sfaira create-dataloader``).
* 1b: Completion of created files (manual).


1a. Create template files.
    .. code-block::

        sfaira create-dataloader --path-data [--path-loader]
    ..
    When creating a dataloader with ``sfaira create-dataloader`` dataloader specific attributes such as organ, organism
    and many more are prompted for.
    We provide a description of all meta data items at the bottom of this page,
    note that these metadata underly specific formattig and ontology constraints described below.
    If the requested information is not available simply hit enter to skip the entry.
    If `--path-loader` is not provided the following default location will be used: `./sfaira/data/dataloaders/loaders/`,
    which is the correct location if you want to commit and push changes from this sfaira clone.
    The CLI decides on an `ID` of this dataset within the loader that you are writing, this will be used to label
    all files associated with the current dataset.
    The CLI tells you how to continue from here, phase 1b) is always necessary, phase 2) is case-dependent and mistakes
    in naming the data folder in phase Pd) are flagged here.
1b. Manual completion of created files (manual).
    1. Correct yaml file.
        Correct errors in `<path_loader>/<DOI-name>/ID.yaml` file and add
        further attributes you may have forgotten in step 2.
        See :ref:`sec-multiple-files` for short-cuts if you have multiple data sets.
        This step is can be skipped if there are the `.yaml` is complete after phase 1a).
    2. Write load function.
        Complete the `load()` function in `<path_loader>/<DOI-name>/ID.py`.

Phase 2: annotate
~~~~~~~~~~~~~~~~~~~

Phase 2 is sub-structured into 2 sub-phases:

* 2a: Create metadata annotation files (``sfaira annotate-dataloader``).
* 2b: Completion of annotation (manual).

Phase 2 can be entirely skipped if no annotation maps are necessary, this is indicated by the CLI at the end of phase 1a.

2a. Create metadata annotation files (``sfaira annotate-dataloader``).
    .. code-block::

        sfaira annotate-dataloader --doi --path_data [--path_loader]
    ..
    ``sfaira annotate-dataloader`` will run fuzzy string matching between the annotations in the metadata column you provided in the
    `cell_types_original_obs_key` attribute of the yaml file and the Cell Ontology Database.
    Note that this will abort with error if there are bugs in your data loader.
2b. Completion of annotation (manual).
    Sfaira creates suggestions for ontology maps in `<path_loader>/<DOI-name>/ID*.tsv` files.
    One such file is created for each meta data item that is annotated per cell,
    as an `ITEM_obs_key` rather than `ITEM` in the `.yaml`.
    This means that a variable number of such files is created and dependending on the scenario, even no such files may
    be necessary.
    Each file contains two columns with one row for each unique free-text meta data item, e.g. each cell type label.

    - The first column is labeled "source" and contains free-text identifiers.
    - The second column is labeled "target" and contains suggestions for matching the symbols from the corresponding ontology.

    The suggestions are based on multiple search criteria, mostly on similarity of the free-text token to tokes in the
    ontology.
    Suggested tokens are separated by ":" in the target column,
    for each token, the same number of suggestions is supplied.
    We use different search strategies on each token and separate the output by strategy by ":||:".
    You might notice that one strategy works well for a particular `ID*.tsv` and focus your attention on that group.
    It is now up to you to manually mitigate the suggestions in the "target" column of each `.tsv` file,
    for example in a text editor.
    Depending on the ontology and on the accuracy of the free-text annotation, these suggestions may be more or
    less helpful.
    The worst case is that you need to go to search engine of the ontology at hand for each entry to check for matches.
    The best case is that you know the ontology well enough to choose from the suggestions,
    assuming that the best match is in the suggestions.
    Reality lies somewhere in the middle of the two, do not be too conservative with looking items up online.
    We suggest to use the ontology search engine on the OLS_ web-interface for your manual queries.

    Note 1: If you compare these `ID*.tsv` to `tsv` files from published data loaders,
    you will notice that published ones contain a third column.
    This column is automatically added in phase 3 if the second column was correctly filled here.

    Note 2: The two columns in the `ID*.tsv` are separated by a tab-separator ("\\t"),
    make sure to not accidentally delete this token.
    If you accidentally replace it with `" "`, you will receive errors in phase 3, so do a visual check after finishing
    your work on each `ID*.tsv` file.

.. _OLS:https://www.ebi.ac.uk/ols/ontologies/cl

Phase 3: finalize
~~~~~~~~~~~~~~~~~~~~

3a. Clean and test data loader.
    .. code-block::

        sfaira finalize-dataloader --doi --path_data [--path_loader]
    ..
    This command will test data loading and will format the metadata maps in `ID*.tsv` files from phase 2b).
    If this command passes without further change requests, the data loader is finished and ready for publication!

Phase 4: upload
~~~~~~~~~~~~~~~~~

4a. Push data loader to the public sfaira repository.
    You can contribute the data loader to public sfaira as code through a pull request.
    Note that you can also just keep the data loader in your local installation if you do not want to make it public.

    .. code-block::

        cd DIR_SFAIRA
        cd sfaira
        git add *
        git commit -m "Completed data loader."
        git push TODO put full line for branch creation here.
    ..

Phase 5: export-h5ad
~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5 and 6 are optional, see also introduction paragraphs on this documentation page.

5a. Export `.h5ads`'s.
    Write streamlined dataset(s) corresponding to data loader into (an) `.h5ad` file(s) according to a specific set of
    rules (a schema).
    .. code-block::

        sfaira export-h5ad --doi --schema --path-out --path_data [--path_loader]
    ..

Phase 6: validate-h5ad
~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5 and 6 are optional, see also introduction paragraphs on this documentation page.

6a. Validate format of `.h5ad` according to a specific set of rules (a schema).
    .. code-block::

        sfaira validate-h5ad --h5ad --schema
    ..



Advanced topics
----------------

.. _sec-multiple-files:
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

Metadata conventions and ontologies
------------------------------------

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
