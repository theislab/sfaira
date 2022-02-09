import os
import pathlib

from sfaira.data import DatasetSuperGroup
from sfaira.data.dataloaders.databases import DatasetSuperGroupDatabases
from sfaira.unit_tests.data_for_tests.databases.consts import CELLXGENE_COLLECTION_ID

from sfaira.unit_tests.directories import DIR_DATA_DATABASES_CACHE


def prepare_dsg_database(database: str, download: bool = True) -> DatasetSuperGroup:
    """
    Prepares data set super group of data base returns instance.

    :param database: Database to make available.
    :param download: Whether to make sure that raw files are downloaded.
    """
    if not os.path.exists(DIR_DATA_DATABASES_CACHE):
        pathlib.Path(DIR_DATA_DATABASES_CACHE).mkdir(parents=True, exist_ok=True)
    if database == "cellxgene":
        # Note: if unit tests related to cellxgene throw error because of caching, empty cache.
        # Unit tests are slow if cache_metadata is False.
        dsg = DatasetSuperGroupDatabases(data_path=DIR_DATA_DATABASES_CACHE, cache_metadata=True)
        # Only retain pre-defined target collections to avoid bulk downloads during unit tests.
        dsg.subset(key="collection_id", values=CELLXGENE_COLLECTION_ID)
    else:
        assert False, database
    if download:
        dsg.download()
    return dsg
