import abc
import anndata
import numpy as np
from typing import Union

from sfaira.estimators import EstimatorKeras, EstimatorKerasCelltype, EstimatorKerasEmbedding
from sfaira.versions.topologies import TopologyContainer

GENES = ["ENSMUSG00000000003", "ENSMUSG00000000028"]
TARGETS = ["T cell", "stromal cell"]

TOPOLOGY_EMBEDDING_MODEL = {
    "model_type": None,
    "input": {
        "genome": "Mus_musculus.GRCm38.102",
        "genes": ["ensg", GENES],
    },
    "output": {},
    "hyper_parameters": {
        "latent_dim": None,
        "l2_coef": 0.,
        "l1_coef": 0.,
        "output_layer": "nb_const_disp"
    }
}

TOPOLOGY_CELLTYPE_MODEL = {
    "model_type": None,
    "input": {
        "genome": "Mus_musculus.GRCm38.102",
        "genes": ["ensg", GENES],
    },
    "output": {
        "cl": "v2021-02-01",
        "targets": TARGETS
    },
    "hyper_parameters": {
        "latent_dim": None,
        "l1_coef": 0.,
        "l2_coef": 0.,
    }
}


class TestEstimatorBase:
    estimator: Union[EstimatorKeras]
    data: Union[anndata.AnnData]

    """
    Contains functions _test* to test individual functions and attributes of estimator class.

    TODO for everybody working on this, add one _test* function in here and add it into
    basic_estimator_test(). See _test_call() for an example.
    """

    def simulate(self):
        """
        Simulate basic data example used for unit test.

        Sets attribute .data with simulated data.

        :return:
        """
        nobs = 100
        self.data = anndata.AnnData(
            np.random.randint(low=0, high=100, size=(nobs, len(GENES))).astype(np.float32)
        )
        self.data.obs["cell_ontology_class"] = [
            TARGETS[np.random.randint(0, len(TARGETS))]
            for i in range(nobs)
        ]
        self.data.var["ensembl"] = GENES

    @abc.abstractmethod
    def init_estimator(self, model_type: str):
        """
        Initialise target estimator as .estimator attribute.
        """
        pass

    @abc.abstractmethod
    def basic_estimator_test(self):
        pass

    def test_for_fatal(self, model_type):
        np.random.seed(1)
        self.simulate()
        self.init_estimator(model_type=model_type)
        self.basic_estimator_test()
        return True


class TestEstimatorKerasEmbedding(TestEstimatorBase):

    estimator: EstimatorKerasEmbedding

    def init_estimator(self, model_type):
        topology = TOPOLOGY_EMBEDDING_MODEL.copy()
        topology["model_type"] = model_type
        if model_type == "linear":
            topology["hyper_parameters"]["latent_dim"] = 2
        else:
            topology["hyper_parameters"]["latent_dim"] = (len(GENES), 2, len(GENES))
        self.model_type = model_type
        self.estimator = EstimatorKerasEmbedding(
            data=self.data,
            model_dir=None,
            model_id="testid",
            model_topology=TopologyContainer(topology=topology, topology_id="0.1")
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
        if self.model_type != 'vae':
            assert np.allclose(prediction_output, new_prediction_output, rtol=1e-6, atol=1e-6)
            assert np.allclose(prediction_embed, new_prediction_embed, rtol=1e-6, atol=1e-6)


class TestEstimatorKerasCelltype(TestEstimatorBase):

    estimator: EstimatorKerasCelltype

    def init_estimator(self, model_type: str):
        topology = TOPOLOGY_CELLTYPE_MODEL.copy()
        topology["model_type"] = model_type
        topology["hyper_parameters"]["latent_dim"] = (len(GENES), 2)
        self.model_type = model_type
        self.estimator = EstimatorKerasCelltype(
            data=self.data,
            model_dir=None,
            model_id="testid",
            model_topology=TopologyContainer(topology=topology, topology_id="0.1"),
        )
        self.estimator.celltype_universe.leaves = TARGETS

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


# Test embedding models:


def test_for_fatal_linear():
    test_estim = TestEstimatorKerasEmbedding()
    test_estim.test_for_fatal(model_type="linear")


def test_for_fatal_ae():
    test_estim = TestEstimatorKerasEmbedding()
    test_estim.test_for_fatal(model_type="ae")


def test_for_fatal_vae():
    test_estim = TestEstimatorKerasEmbedding()
    test_estim.test_for_fatal(model_type="vae")


# Test cell type predictor models:


def test_for_fatal_mlp():
    test_estim = TestEstimatorKerasCelltype()
    test_estim.test_for_fatal(model_type="mlp")


def test_for_fatal_marker():
    test_estim = TestEstimatorKerasCelltype()
    test_estim.test_for_fatal(model_type="marker")
