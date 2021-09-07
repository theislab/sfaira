import os
import pytest
import shutil
from typing import List

from sfaira.unit_tests.directories import DIR_DATA_DATABASES_CACHE
from sfaira.unit_tests.data_for_tests.databases.utils import prepare_dsg_database
from sfaira.unit_tests.data_for_tests.databases.consts import CELLXGENE_DATASET_ID


# Execute this one first so that data sets are only downloaded once. Named test_a for this reason.
@pytest.mark.parametrize("database", ["cellxgene", ])
@pytest.mark.parametrize("subset_args", [None, ["id", CELLXGENE_DATASET_ID], ])
def test_a_dsgs_download(database: str, subset_args: List[str]):
    """
    Tests if downloading of data base entries works.

    Warning, deletes entire database unit test cache.
    """
    if os.path.exists(DIR_DATA_DATABASES_CACHE):
        shutil.rmtree(DIR_DATA_DATABASES_CACHE)
    dsg = prepare_dsg_database(database=database, download=False)
    if subset_args is not None:
        dsg.subset(key=subset_args[0], values=subset_args[1])
    dsg.download()


@pytest.mark.parametrize("database", ["cellxgene", ])
@pytest.mark.parametrize("subset_args", [["id", CELLXGENE_DATASET_ID], ["organism", "human"], ])
def test_dsgs_subset(database: str, subset_args: List[str]):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    dsg = prepare_dsg_database(database=database)
    dsg.subset(key=subset_args[0], values=subset_args[1])
