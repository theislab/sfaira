import abc
import numpy as np
import os
import pandas as pd
from typing import Union
import unittest

from sfaira.interface.model_zoo import ModelZoo, ModelZooCelltype, ModelZooEmbedding


class _TestZoo:
    zoo: Union[ModelZoo]
    data: np.ndarray

    """
    Contains functions _test* to test individual functions and attributes of estimator class.

    TODO for everybody working on this, add one _test* function in here and add it into
    basic_estimator_test(). See _test_kipoi_call() for an example.
    """

    @abc.abstractmethod
    def init_zoo(self):
        """
        Initialise target zoo as .zoo attribute.

        :return:
        """
        pass

    def simulate(self):
        """
        Simulate basic data example used for unit test.

        Sets attribute .data with simulated data.

        :return:
        """
        pass

    def _test_kipoi_call(self):
        """
        Test whether kipoi_experimental model call works.

        :return:
        """
        self.zoo.call_kipoi()

    def _test_basic(self, id: str):
        """
        Test all relevant model methods.


        :return:
        """
        np.random.seed(1)
        self.simulate()
        self.init_zoo()
        # self._test_kipoi_call()
        self.zoo_manual.set_model_id(id)


class TestZooKerasEmbedding(unittest.TestCase, _TestZoo):

    def init_zoo(self):
        package_dir = str(os.path.dirname(os.path.abspath(__file__)))
        lookup_table = pd.read_csv(
            os.path.join(package_dir, '../test_data', 'model_lookuptable.csv'),
            header=0, index_col=0
        )
        self.zoo = ModelZooEmbedding(model_lookuptable=lookup_table)
        self.zoo_manual = ModelZooEmbedding(model_lookuptable=None)

    def test_basic(self):
        self._test_basic(id="embedding_mouse_lung_vae_theislab_0.1_0.1")
        self.zoo.set_latest('mouse', 'lung', 'vae', 'theislab', '0.1')
        assert self.zoo.model_id == "embedding_mouse_lung_vae_theislab_0.1_0.1"
        assert self.zoo.model_id == self.zoo_manual.model_id


class TestZooKerasCelltype(unittest.TestCase, _TestZoo):

    def init_zoo(self):
        package_dir = str(os.path.dirname(os.path.abspath(__file__)))
        lookup_table = pd.read_csv(
            os.path.join(package_dir, '../test_data', 'model_lookuptable.csv'),
            header=0, index_col=0
        )
        self.zoo = ModelZooCelltype(model_lookuptable=lookup_table)
        self.zoo_manual = ModelZooCelltype(model_lookuptable=None)

    def test_basic(self):
        self._test_basic(id="celltype_mouse_lung_mlp_theislab_0.0.1_0.1")
        self.zoo.set_latest('mouse', 'lung', 'mlp', 'theislab', '0.0.1')
        assert self.zoo.model_id == "celltype_mouse_lung_mlp_theislab_0.0.1_0.1"
        assert self.zoo.model_id == self.zoo_manual.model_id


if __name__ == '__main__':
    unittest.main()
