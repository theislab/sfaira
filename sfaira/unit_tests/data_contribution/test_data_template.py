import os
import pydoc
import shutil

from sfaira.data import DatasetGroupDirectoryOriented

try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None


def _get_ds(doi_sfaira_repr: str, test_data: str):
    dir_loader_sfaira = "sfaira.data.dataloaders.loaders."
    file_path_sfaira = "/" + "/".join(pydoc.locate(dir_loader_sfaira + "FILE_PATH").split("/")[:-1])
    if sfairae is not None:
        dir_loader_sfairae = "sfaira_extension.data.dataloaders.loaders."
        file_path_sfairae = "/" + "/".join(pydoc.locate(dir_loader_sfairae + "FILE_PATH").split("/")[:-1])
    else:
        file_path_sfairae = None
    # Check if loader name is a directory either in sfaira or sfaira_extension loader collections:
    if doi_sfaira_repr in os.listdir(file_path_sfaira):
        dir_loader = dir_loader_sfaira + "." + doi_sfaira_repr
    elif doi_sfaira_repr in os.listdir(file_path_sfairae):
        dir_loader = dir_loader_sfairae + "." + doi_sfaira_repr
    else:
        raise ValueError("data loader not found in sfaira and also not in sfaira_extension")
    file_path = pydoc.locate(dir_loader + ".FILE_PATH")
    # Clear dataset cache
    cache_path = test_data
    ds = DatasetGroupDirectoryOriented(
        file_base=file_path,
        data_path=test_data,
        meta_path=None,
        cache_path=cache_path
    )
    return ds, cache_path


def test_load(doi_sfaira_repr: str, test_data: str):
    ds, cache_path = _get_ds(doi_sfaira_repr=doi_sfaira_repr, test_data=test_data)
    ds.clean_ontology_class_map()
    ds.load(
        remove_gene_version=True,
        load_raw=True,
        allow_caching=True,
        verbose=3,
    )
    # Try loading from cache:
    ds, _ = _get_ds(doi_sfaira_repr=doi_sfaira_repr, test_data=test_data)
    # TODO try-except with good error description saying that the data loader is broken here:
    ds.load(
        remove_gene_version=True,
        load_raw=False,
        allow_caching=True,
        verbose=3,
    )
    shutil.rmtree(os.path.join(cache_path, doi_sfaira_repr, "cache"), ignore_errors=True)
