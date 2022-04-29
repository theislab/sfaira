import os
import shutil
from typing import Union

from sfaira import settings


def clean_cache(cache: Union[None, str] = None):
    """
    Utility function to clean cached objects in paths of sfaira installation.

    This can be used to force re-caching or to reduce directory size.
    """
    if cache is not None:
        cache_dir_dict = {
            "all": settings.cachedir_base,
            "dataset_meta": settings.cachedir_databases,
            "genomes": settings.cachedir_genomes,
            "ontologies": settings.cachedir_ontologies,
        }
        if cache not in cache_dir_dict.keys():
            raise ValueError(f"Did not find cache directory input {cache} in support list: "
                             f"{list(cache_dir_dict.keys())}")
        else:
            print(f"cleaning cache {cache} in directory {cache_dir_dict[cache]}")
            # Assert that a path within sfaira is selected as a sanity check:
            dir_to_delete = cache_dir_dict[cache]
            dir_sfaira = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            assert str(dir_to_delete).startswith(dir_sfaira), \
                f"trying to delete outside of sfaira installation: {dir_to_delete}"
            shutil.rmtree(dir_to_delete)


def clean_doi(doi: str):
    # Check if DOI is already in sfaira format, in that case do not prefix with "d":
    if doi.startswith("d10_") or doi.startswith("dno_doi_"):
        return f'{doi.translate({ord(c): "_" for c in r"!@#$%^&*()[]/{};:,.<>?|`~-=_+"})}'
    else:
        return f'd{doi.translate({ord(c): "_" for c in r"!@#$%^&*()[]/{};:,.<>?|`~-=_+"})}'


def clean_id_str(s):
    if s is not None:
        s = s.replace(',', '').replace(' ', '').replace('-', '').replace('_', '').replace("'", '').lower()
    return s
