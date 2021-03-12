import os
import pydoc

from sfaira.data import DatasetGroupDirectoryOriented, DatasetGroup, DatasetBase
from sfaira.data.utils import read_yaml
try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None


def test_load(dir_template: str = "../template_data", doi_sfaira_repr="d10_1016_j_cmet_2019_01_021"):
    """
    Unit test to assist with data set contribution.

    The workflow for contributing a data set with this data loader is as follows:

    1. Write a data loader and add it into the loader directory of your local sfaira installation.
    2. Address ToDos below.
    3. Run this unit test until you are not getting errors from your data loader anymore.

    In the process of this unit test, this data loader will have written putative cell type maps from your
    annotation to the cell ontology.

    4. Moderate the suggestions made here: Choose the best fit cell ontology label for your cells.
    Sfaira uses multiple mechanisms of finding matches, depending on how the free text was generated, these might be
    differentially successful. The proposed IDs groups are separate by ":|||:" strings to give you a visual anchor
    when going through these lists. You need to delete all of these division strings and all labels in the second
    columns other than the best fit label. Do not change the first column,
    (Note that columns are separated by ",")
    You can also manually check maps here: https://www.ebi.ac.uk/ols/ontologies/cl
    5. Run this unit test for a last time to check the cell type maps.

    :return:
    """
    remove_gene_version = True
    match_to_reference = None

    flattened_doi = doi_sfaira_repr
    # Define file names and loader paths in sfaira or sfaira_extension:
    # Define base paths of loader collections in sfaira and sfaira_extension:
    dir_loader_sfaira = "sfaira.data.dataloaders.loaders."
    file_path_sfaira = "/" + "/".join(pydoc.locate(dir_loader_sfaira + "FILE_PATH").split("/")[:-1])
    if sfairae is not None:
        dir_loader_sfairae = "sfaira_extension.data.dataloaders.loaders."
        file_path_sfairae = "/" + "/".join(pydoc.locate(dir_loader_sfairae + "FILE_PATH").split("/")[:-1])
    else:
        file_path_sfairae = None
    # Check if loader name is a directory either in sfaira or sfaira_extension loader collections:
    if flattened_doi in os.listdir(file_path_sfaira):
        dir_loader = dir_loader_sfaira + "." + flattened_doi
        package_source = "sfaira"
    elif flattened_doi in os.listdir(file_path_sfairae):
        dir_loader = dir_loader_sfairae + "." + flattened_doi
        package_source = "sfairae"
    else:
        raise ValueError("data loader not found in sfaira and also not in sfaira_extension")
    file_path = pydoc.locate(dir_loader + ".FILE_PATH")

    ds = DatasetGroupDirectoryOriented(
        file_base=file_path,
        data_path=dir_template,
        meta_path=dir_template,
        cache_path=dir_template
    )
    # Test raw loading and caching:
    # You can set load_raw to True while debugging when caching works already to speed the test up,
    # but be sure to set load_raw to True for final tests.
    ds.load(
        remove_gene_version=False,
        match_to_reference=False,
        load_raw=True,  # tests raw loading
        allow_caching=True,  # tests caching
        set_metadata=False,
    )
    assert len(ds.ids) > 0, f"no data sets loaded, make sure raw data is in {dir_template}"
    # Create cell type conversion table:
    cwd = os.path.dirname(file_path)
    dataset_module = str(cwd.split("/")[-1])
    # Group data sets by file module:
    # Note that if we were not grouping the cell type map .tsv files by file module, we could directly call
    # write_ontology_class_map on the ds.
    for f in os.listdir(cwd):
        if os.path.isfile(os.path.join(cwd, f)):  # only files
            # Narrow down to data set files:
            if f.split(".")[-1] == "py" and f.split(".")[0] not in ["__init__", "base", "group"]:
                file_module = ".".join(f.split(".")[:-1])

                # I) Instantiate Data set group to get all IDs of data sets associated with this .py file.
                # Note that all data sets in this directory are already loaded in ds, so we just need the IDs.
                DatasetFound = pydoc.locate(dir_loader + "." + file_module + ".Dataset")
                # Load objects from name space:
                # - load(): Loading function that return anndata instance.
                # - SAMPLE_FNS: File name list for DatasetBaseGroupLoadingManyFiles
                load_func = pydoc.locate(dir_loader + "." + file_module + ".load")
                load_func_annotation = pydoc.locate(dir_loader + "." + file_module + ".LOAD_ANNOTATION")
                # Also check sfaira_extension for additional load_func_annotation:
                if package_source != "sfairae":
                    load_func_annotation_sfairae = pydoc.locate(dir_loader_sfairae + "." + dataset_module +
                                                                "." + file_module + ".LOAD_ANNOTATION")
                    # LOAD_ANNOTATION is a dictionary so we can use update to extend it.
                    if load_func_annotation_sfairae is not None and load_func_annotation is not None:
                        load_func_annotation.update(load_func_annotation_sfairae)
                    elif load_func_annotation_sfairae is not None and load_func_annotation is None:
                        load_func_annotation = load_func_annotation_sfairae
                sample_fns = pydoc.locate(dir_loader + "." + file_module + ".SAMPLE_FNS")
                fn_yaml = os.path.join(cwd, file_module + ".yaml")
                fn_yaml = fn_yaml if os.path.exists(fn_yaml) else None
                # Check for sample_fns in yaml:
                if fn_yaml is not None:
                    assert os.path.exists(fn_yaml), f"did not find yaml {fn_yaml}"
                    yaml_vals = read_yaml(fn=fn_yaml)
                    if sample_fns is None and yaml_vals["meta"]["sample_fns"] is not None:
                        sample_fns = yaml_vals["meta"]["sample_fns"]
                if sample_fns is None:
                    sample_fns = [None]
                # Here we distinguish between class that are already defined and those that are not.
                # The latter case arises if meta data are defined in YAMLs and _load is given as a function.
                if DatasetFound is None:
                    datasets_f = [
                        DatasetBase(
                            data_path=dir_template,
                            meta_path=dir_template,
                            cache_path=dir_template,
                            load_func=load_func,
                            dict_load_func_annotation=load_func_annotation,
                            sample_fn=x,
                            sample_fns=sample_fns if sample_fns != [None] else None,
                            yaml_path=fn_yaml,
                        ) for x in sample_fns
                    ]
                else:
                    datasets_f = [
                        DatasetFound(
                            data_path=dir_template,
                            meta_path=dir_template,
                            cache_path=dir_template,
                            load_func=load_func,
                            load_func_annotation=load_func_annotation,
                            sample_fn=x,
                            sample_fns=sample_fns if sample_fns != [None] else None,
                            yaml_path=fn_yaml,
                        ) for x in sample_fns
                    ]
                # II) Build a data set group from the already loaded data sets and use the group ontology writing
                # function.
                dsg_f = DatasetGroup(datasets=dict([(x.id, ds.datasets[x.id]) for x in datasets_f]))
                # III) Write this directly into sfaira installation so that it can be committed via git.
                dsg_f.write_ontology_class_map(
                    fn=os.path.join(cwd, file_module + ".tsv"),
                    protected_writing=True,
                    n_suggest=4,
                )

    # Test loading from cache:
    ds = DatasetGroupDirectoryOriented(
        file_base=file_path,
        data_path=dir_template,
        meta_path=dir_template,
        cache_path=dir_template
    )
    ds.load(
        remove_gene_version=remove_gene_version,
        match_to_reference=match_to_reference,
        load_raw=False,
        allow_caching=False
    )
    ds.clean_ontology_class_map()
    # Test concatenation:
    _ = ds.adata
