import numpy as np
import os
from typing import Union
import pandas as pd
import urllib.request

from sfaira import settings
from sfaira.ui import UserInterface
from sfaira.unit_tests.data_for_tests.loaders.utils import PrepareData
from sfaira.unit_tests import DIR_TEMP


class HelperUi:
    ui: Union[UserInterface]
    data: np.ndarray

    """
    Contains functions _test* to test individual functions and attributes of the user ui class.

    TODO for everybody working on this, add one _test* function in here and add it into
    basic_estimator_test(). See _test_call() for an example.
    """

    def __init__(self):
        self.temp_fn = os.path.join(DIR_TEMP, "test_data")

    def prepare_local_tempfiles(self):
        # create temp directory
        if not os.path.exists(self.temp_fn):
            os.makedirs(self.temp_fn)
        # download an example weight from sfaira repo
        lookuptable = pd.read_csv(
            os.path.join(settings.sfaira_repo_url, 'model_lookuptable.csv'),
            header=0,
            index_col=0
        )
        url = lookuptable.loc[0, "model_file_path"]
        if os.path.basename(url) not in os.listdir(self.temp_fn):
            urllib.request.urlretrieve(url, os.path.join(self.temp_fn, os.path.basename(url)))

    def _get_adata(self):
        """
        Create an adata object for use in unit tests

        :return:
        """
        dsg = PrepareData().prepare_dsg(rewrite=True, load=False)
        dsg.subset(key="id", values=["homosapiens_lung_2021_None_mock4_001_no_doi_mock4"])
        dsg.load()
        return dsg.adata

    def test_local_repo_ui_init(self):
        """
        Test if the sfaira UI can be sucessfully initialised using a local model repository

        :return:
        """
        self.ui = UserInterface(custom_repo=self.temp_fn, sfaira_repo=False)

    def test_public_repo_ui_init(self):
        """
        Test if the sfaira UI can be sucessfully initialised using the public sfaira model repository

        :return:
        """
        self.ui = UserInterface(custom_repo=None, sfaira_repo=True, cache_path=self.temp_fn)

    def test_data_and_model_loading(self):
        self.ui = UserInterface(custom_repo=None, sfaira_repo=True, cache_path=self.temp_fn)
        self.ui.zoo_embedding.model_id = 'embedding_human-blood-ae-0.2-0.1_theislab'
        self.ui.zoo_celltype.model_id = 'celltype_human-blood-mlp-0.1.3-0.1_theislab'
        test_data = self._get_adata()
        self.ui.load_data(test_data, gene_ens_col='index')
        self.ui.load_model_celltype()
        self.ui.load_model_embedding()
        self.ui.predict_all()
        assert "X_sfaira" in self.ui.data.adata.obsm_keys()
        assert "celltypes_sfaira" in self.ui.data.adata.obs_keys()


def test_ui():
    ui = HelperUi()
    ui.prepare_local_tempfiles()
    ui.test_public_repo_ui_init()
    ui.test_local_repo_ui_init()
    ui.test_data_and_model_loading()
