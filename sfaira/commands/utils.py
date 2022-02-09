import os
import re
import pydoc
import sys

from sfaira.commands.consts import PACKAGE_LOADER_PATH, PACKAGE_SFAIRAE_LOADER_PATH

try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None


def doi_lint(x: str) -> bool:
    """
    Check if string matches DOI format.

    Allowed formats are:
        - an actual DOI startin with "10*"
        - a string starting with "no_doi*" indicating that a dataset without corresponding DOI is used
        - a sfaira formated DOI starting with "d10*" or "dno_doi*"
    """
    return re.match(r'\b10\.\d+/[\w.]+\b', x) or x.startswith("no_doi") or \
        x.startswith("d10_") or x.startswith("dno_doi")


def get_pydoc(path_loader, doi_sfaira_repr) -> [str, str]:
    pydoc_handle_sfaira = "sfaira.data.dataloaders.loaders"
    file_path_sfaira = os.path.dirname(str(pydoc.locate(pydoc_handle_sfaira + ".FILE_PATH")))
    pydoc_handle_sfairae = "sfaira_extension.data.dataloaders.loaders" if sfairae else None
    file_path_sfairae = os.path.dirname(str(pydoc.locate(pydoc_handle_sfairae + ".FILE_PATH"))) if sfairae else None

    # Check if loader name is a directory either in external, sfaira or sfaira_extension loader collections:
    if path_loader != PACKAGE_LOADER_PATH and path_loader != PACKAGE_SFAIRAE_LOADER_PATH:  # external
        if doi_sfaira_repr in os.listdir(path_loader):
            if path_loader not in sys.path:
                sys.path.append(path_loader)
            pydoc_handle = doi_sfaira_repr
            # Imitate result from pydoc locate:
            file_path = os.path.join(path_loader, doi_sfaira_repr, "__init__.py")
        else:
            raise ValueError(f"did not find loader {doi_sfaira_repr} in custom path {os.listdir(path_loader)}")
    elif doi_sfaira_repr in os.listdir(file_path_sfaira):
        pydoc_handle = pydoc_handle_sfaira + "." + doi_sfaira_repr
        file_path = os.path.join(file_path_sfaira, doi_sfaira_repr, "__init__.py")
    elif file_path_sfairae and doi_sfaira_repr in os.listdir(file_path_sfairae):
        pydoc_handle = pydoc_handle_sfairae + "." + doi_sfaira_repr
        file_path = os.path.join(file_path_sfairae, doi_sfaira_repr, "__init__.py")
    else:
        raise ValueError("data loader not found in sfaira and also not in sfaira_extension")
    return file_path, pydoc_handle
