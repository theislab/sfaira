"""
All paths used throughout unit testing for temporary files.
"""

import os
import shutil

DIR_BASE_TESTS = os.path.dirname(__file__)
DIR_MODULE_ROOT = os.path.dirname(DIR_BASE_TESTS)
DIR_TEMP = os.path.join(DIR_BASE_TESTS, "temp")

_DIR_DATA_LOADERS = os.path.join(DIR_TEMP, "loaders")
DIR_DATA_LOADERS_CACHE = os.path.join(_DIR_DATA_LOADERS, "cache")
DIR_DATA_LOADERS_STORE_DAO = os.path.join(_DIR_DATA_LOADERS, "store_dao")
DIR_DATA_LOADERS_STORE_H5AD = os.path.join(_DIR_DATA_LOADERS, "store_h5ad")
_DIR_DATA_DATABASES = os.path.join(DIR_TEMP, "databases")
DIR_DATA_DATABASES_CACHE = os.path.join(_DIR_DATA_DATABASES, "cache")
DIR_DATABASE_STORE_DAO = os.path.join(_DIR_DATA_DATABASES, "store_dao")

DIR_SFAIRA_LOADERS = os.path.join(DIR_MODULE_ROOT, "sfaira", "data", "dataloaders", "loaders")


def save_delete(fn):
    assert str(fn).startswith(DIR_TEMP), f"tried to delete outside of temp directory {fn}"
    if os.path.isdir(fn):
        shutil.rmtree(fn)
    else:
        os.remove(fn)
