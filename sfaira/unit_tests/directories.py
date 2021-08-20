"""
All paths used throughout unit testing for temporary files.
"""

import os

DIR_TEMP = os.path.join(os.path.dirname(__file__), "temp")

_DIR_DATA_LOADERS = os.path.join(DIR_TEMP, "loaders")
DIR_DATA_LOADERS_CACHE = os.path.join(_DIR_DATA_LOADERS, "cache")
DIR_DATA_LOADERS_STORE_DAO = os.path.join(_DIR_DATA_LOADERS, "store_dao")
DIR_DATA_LOADERS_STORE_H5AD = os.path.join(_DIR_DATA_LOADERS, "store_h5ad")
_DIR_DATA_DATABASES = os.path.join(DIR_TEMP, "databases")
DIR_DATA_DATABASES_CACHE = os.path.join(_DIR_DATA_DATABASES, "cache")
