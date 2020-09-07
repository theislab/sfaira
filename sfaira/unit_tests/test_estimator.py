import abc
import anndata
import numpy as np
import tensorflow as tf
from typing import Union
import unittest

from sfaira.unit_tests.external import EstimatorKeras, EstimatorKerasCelltype, EstimatorKerasEmbedding
from sfaira.unit_tests.external import celltype_versions, SuperGenomeContainer, Topologies


class _TestEstimator:
    estimator: Union[EstimatorKeras]
    data: Union[anndata.AnnData]
    model: Union[tf.keras.models.Model, None]
    topology_container: Union[Topologies, None]
    optimizer: Union[str, None]
    model_id: Union[str, None]
    weights: Union[np.ndarray, None]
    model_dir: Union[str, None]

    """
    Contains functions _test* to test individual functions and attributes of estimator class.
    
    TODO for everybody working on this, add one _test* function in here and add it into 
    basic_estimator_test(). See _test_call() for an example.
    """

    @abc.abstractmethod
    def set_topology(self, model_type):
        pass

    def simulate(self):
        """
        Simulate basic data example used for unit test.

        Sets attribute .data with simulated data.

        :return:
        """
        nobs = 100
        ngenes = self.topology_container.ngenes
        self.data = anndata.AnnData(
            np.random.randint(low=0, high=100, size=(nobs, ngenes)).astype(np.float32)
        )
        self.data.obs["cell_ontology_class"] = [
            ["vein endothelial cell", "glial cell"][np.random.randint(0, 2)]
            for i in range(nobs)
        ]
        self.data.var["ensembl"] = self.topology_container.genome_container.ensembl
        self.data.var["names"] = self.topology_container.genome_container.names

    @abc.abstractmethod
    def init_estimator(self):
        """
        Initialise target estimator as .estimator attribute.
        """
        pass

    @abc.abstractmethod
    def basic_estimator_test(self):
        pass

    def _test_for_fatal(self):
        np.random.seed(1)
        self.simulate()
        self.init_estimator()
        self.basic_estimator_test()
        return True


class TestEstimatorKerasEmbedding(unittest.TestCase, _TestEstimator):

    def set_topology(self, model_type):
        self.topology_container = Topologies(
            species="mouse",
            model_class="embedding",
            model_type=model_type,
            topology_id="0.1"
        )

    def init_estimator(self):
        self.estimator = EstimatorKerasEmbedding(
            data=self.data,
            model_dir=None,
            model_id=None,
            species="mouse",
            organ="lung",
            model_type=self.topology_container.model_type,
            model_topology=self.topology_container.topology_id
        )

    def basic_estimator_test(self):
        self.estimator.init_model()
        self.estimator.train(
            optimizer="adam",
            lr=0.005,
            epochs=2,
            batch_size=32,
            validation_split=0.1,
            test_split=0.1,
            validation_batch_size=32,
            max_validation_steps=1
        )
        _ = self.estimator.evaluate()
        prediction_output = self.estimator.predict()
        prediction_embed = self.estimator.predict_embedding()
        weights = self.estimator.model.training_model.get_weights()
        self.estimator.save_weights_to_cache()
        self.estimator.load_weights_from_cache()
        new_prediction_output = self.estimator.predict()
        new_prediction_embed = self.estimator.predict_embedding()
        new_weights = self.estimator.model.training_model.get_weights()
        for i in range(len(weights)):
            assert np.allclose(weights[i], new_weights[i], rtol=1e-6, atol=1e-6)
        if self.topology_container.model_type != 'vae':
            assert np.allclose(prediction_output, new_prediction_output, rtol=1e-6, atol=1e-6)
            assert np.allclose(prediction_embed, new_prediction_embed, rtol=1e-6, atol=1e-6)

    def test_for_fatal_vae(self):
        self.set_topology(model_type="vae")
        self._test_for_fatal()

    def test_for_fatal_ae(self):
        self.set_topology(model_type="ae")
        self._test_for_fatal()

    def test_for_fatal_linear(self):
        self.set_topology(model_type="linear")
        self._test_for_fatal()


class TestEstimatorKerasCelltype(unittest.TestCase, _TestEstimator):

    def set_topology(self, model_type):
        self.topology_container = Topologies(
            species="mouse",
            model_class="celltype",
            model_type=model_type,
            topology_id="0.0.1"
        )

    def init_estimator(self):
        self.estimator = EstimatorKerasCelltype(
            data=self.data,
            model_dir=None,
            model_id=None,
            species="mouse",
            organ="lung",
            model_type=self.topology_container.model_type,
            model_topology=self.topology_container.topology_id
        )

    def basic_estimator_test(self):
        self.estimator.init_model()
        self.estimator.train(
            optimizer="adam",
            lr=0.005,
            epochs=2,
            batch_size=32,
            validation_split=0.1,
            test_split=0.1,
            validation_batch_size=32,
            max_validation_steps=1
        )
        _ = self.estimator.evaluate()
        prediction_output = self.estimator.predict()
        weights = self.estimator.model.training_model.get_weights()
        self.estimator.save_weights_to_cache()
        self.estimator.load_weights_from_cache()
        new_prediction_output = self.estimator.predict()
        new_weights = self.estimator.model.training_model.get_weights()
        print(self.estimator.model.training_model.summary())
        for i in range(len(weights)):
            assert np.allclose(weights[i], new_weights[i], rtol=1e-6, atol=1e-6)
        assert np.allclose(prediction_output, new_prediction_output, rtol=1e-6, atol=1e-6)

    def test_for_fatal_mlp(self):
        self.set_topology(model_type="mlp")
        self._test_for_fatal()

    def test_for_fatal_marker(self):
        self.set_topology(model_type="marker")
        self._test_for_fatal()


if __name__ == '__main__':
    unittest.main()
