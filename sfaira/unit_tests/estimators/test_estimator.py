import abc
import anndata
import numpy as np
import os
import pandas as pd
import pytest
import time
from typing import Union

from sfaira.data import DistributedStore
from sfaira.estimators import EstimatorKeras, EstimatorKerasCelltype, EstimatorKerasEmbedding
from sfaira.versions.topologies import TopologyContainer
from sfaira.unit_tests.utils import cached_store_writing, simulate_anndata

dir_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data")
dir_meta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data/meta")
cache_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
                         "cache", "genomes")

ASSEMBLY = "Mus_musculus.GRCm38.102"
GENES = ["ENSMUSG00000000003", "ENSMUSG00000000028"]
TARGETS = ["T cell", "stromal cell"]
ASSAYS = ["10x sequencing", "Smart-seq2"]


TOPOLOGY_EMBEDDING_MODEL = {
    "model_type": None,
    "input": {
        "genome": ASSEMBLY,
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
        "genome": ASSEMBLY,
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


class HelperEstimatorBase:

    data: Union[anndata.AnnData, DistributedStore]
    estimator: Union[EstimatorKeras]
    model_type: str
    tc: TopologyContainer

    """
    Contains functions _test* to test individual functions and attributes of estimator class.

    TODO for everybody working on this, add one _test* function in here and add it into
    basic_estimator_test(). See _test_call() for an example.
    """

    def _simulate(self) -> anndata.AnnData:
        """
        Simulate basic data example used for unit test.

        :return: Simulated data set.
        """
        return simulate_anndata(n_obs=100, assays=ASSAYS, genes=self.tc.gc.ensembl, targets=TARGETS)

    def load_adata(self):
        """
        Sets attribute .data with simulated data.
        """
        self.data = self._simulate()

    def load_store(self):
        store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=ASSEMBLY)
        store = DistributedStore(cache_path=store_path)
        self.data = store

    @abc.abstractmethod
    def init_topology(self, model_type: str, feature_space: str):
        pass

    @abc.abstractmethod
    def init_estimator(self):
        """
        Initialise target estimator as .estimator attribute.
        """
        pass

    def estimator_train(self, test_split):
        self.estimator.init_model()
        self.estimator.train(
            optimizer="adam",
            lr=0.005,
            epochs=2,
            batch_size=4,
            validation_split=0.5,
            test_split=test_split,
            validation_batch_size=4,
            max_validation_steps=1,
            shuffle_buffer_size=10,
            cache_full=False,
        )

    @abc.abstractmethod
    def basic_estimator_test(self, test_split):
        pass

    def load_estimator(self, model_type, data_type, feature_space, test_split):
        self.init_topology(model_type=model_type, feature_space=feature_space)
        np.random.seed(1)
        if data_type == "adata":
            self.load_adata()
        else:
            self.load_store()
        self.init_estimator()
        self.estimator_train(test_split=test_split)

    def fatal_estimator_test(self, model_type, data_type, test_split=0.1, feature_space="small"):
        self.load_estimator(model_type=model_type, data_type=data_type, feature_space=feature_space,
                            test_split=test_split)
        self.basic_estimator_test()


class HelperEstimatorKerasEmbedding(HelperEstimatorBase):

    estimator: EstimatorKerasEmbedding
    model_type: str
    tc: TopologyContainer

    def init_topology(self, model_type: str, feature_space: str):
        topology = TOPOLOGY_EMBEDDING_MODEL.copy()
        if feature_space == "full":
            # Read 500 genes (not full protein coding) to compromise between being able to distinguish observations
            # and reducing run time of unit tests.
            tab = pd.read_csv(os.path.join(cache_dir, ASSEMBLY + ".csv"))
            genes_full = tab.loc[tab["gene_biotype"].values == "protein_coding", "gene_id"].values[:500].tolist()
            topology["input"]["genes"] = ["ensg", genes_full]
        topology["model_type"] = model_type
        if model_type == "linear":
            topology["hyper_parameters"]["latent_dim"] = 2
        else:
            topology["hyper_parameters"]["latent_dim"] = (2, 2, 2)
        self.model_type = model_type
        self.tc = TopologyContainer(topology=topology, topology_id="0.1")

    def init_estimator(self):
        self.estimator = EstimatorKerasEmbedding(
            data=self.data,
            model_dir=None,
            model_id="testid",
            model_topology=self.tc
        )

    def basic_estimator_test(self, test_split):
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
            if not np.any(np.isnan(weights[i])):
                assert np.allclose(weights[i], new_weights[i], rtol=1e-6, atol=1e-6)
        if self.model_type != 'vae':
            if not np.any(np.isnan(prediction_output)):
                assert np.allclose(prediction_output, new_prediction_output, rtol=1e-6, atol=1e-6)
                assert np.allclose(prediction_embed, new_prediction_embed, rtol=1e-6, atol=1e-6)


class HelperEstimatorKerasCelltype(HelperEstimatorBase):

    estimator: EstimatorKerasCelltype
    model_type: str
    tc: TopologyContainer

    def init_topology(self, model_type: str, feature_space: str):
        topology = TOPOLOGY_CELLTYPE_MODEL.copy()
        topology["model_type"] = model_type
        topology["hyper_parameters"]["latent_dim"] = (2,)
        self.model_type = model_type
        self.tc = TopologyContainer(topology=topology, topology_id="0.1")

    def init_estimator(self):
        self.estimator = EstimatorKerasCelltype(
            data=self.data,
            model_dir=None,
            model_id="testid",
            model_topology=self.tc
        )
        self.estimator.celltype_universe.leaves = TARGETS

    def basic_estimator_test(self, test_split):
        _ = self.estimator.evaluate()
        prediction_output = self.estimator.predict()
        weights = self.estimator.model.training_model.get_weights()
        self.estimator.save_weights_to_cache()
        self.estimator.load_weights_from_cache()
        new_prediction_output = self.estimator.predict()
        new_weights = self.estimator.model.training_model.get_weights()
        print(self.estimator.model.training_model.summary())
        for i in range(len(weights)):
            if not np.any(np.isnan(weights[i])):
                assert np.allclose(weights[i], new_weights[i], rtol=1e-6, atol=1e-6)
        if not np.any(np.isnan(prediction_output)):
            assert np.allclose(prediction_output, new_prediction_output, rtol=1e-6, atol=1e-6)


# Test embedding models:


@pytest.mark.parametrize("data_type", ["adata", "store"])
def test_for_fatal_linear(data_type):
    test_estim = HelperEstimatorKerasEmbedding()
    test_estim.fatal_estimator_test(model_type="linear", data_type=data_type)


def test_for_fatal_ae():
    test_estim = HelperEstimatorKerasEmbedding()
    test_estim.fatal_estimator_test(model_type="ae", data_type="adata")


def test_for_fatal_vae():
    test_estim = HelperEstimatorKerasEmbedding()
    test_estim.fatal_estimator_test(model_type="vae", data_type="adata")


# Test cell type predictor models:


@pytest.mark.parametrize("data_type", ["adata", "store"])
def test_for_fatal_mlp(data_type):
    test_estim = HelperEstimatorKerasCelltype()
    test_estim.fatal_estimator_test(model_type="mlp", data_type=data_type)


def test_for_fatal_marker():
    test_estim = HelperEstimatorKerasCelltype()
    test_estim.fatal_estimator_test(model_type="marker", data_type="adata")


# Test index sets


@pytest.mark.parametrize("data_type", ["adata", "store"])
@pytest.mark.parametrize("test_split", [0.3, {"assay_sc": "10x sequencing"}])
def test_split_index_sets(data_type: str, test_split):
    """
    Test that train, val, test split index sets are correct:

        1) complete
        2) non-overlapping
        3) that test indices map to all (float split) or distinct (attribute split) data sets
        4) do not contain duplicated observations within and across splits (defined based on the feature vectors)
    """
    test_estim = HelperEstimatorKerasEmbedding()
    # Need full feature space here because observations are not necessarily different in small model testing feature
    # space with only two genes:
    t0 = time.time()
    test_estim.load_estimator(model_type="linear", data_type=data_type, test_split=test_split, feature_space="full")
    print(f"time for running estimator test: {time.time() - t0}s")
    idx_train = test_estim.estimator.idx_train
    idx_eval = test_estim.estimator.idx_eval
    idx_test = test_estim.estimator.idx_test
    # 1) Assert that index assignments sum up to full data set:
    assert len(idx_train) + len(idx_eval) + len(idx_test) == test_estim.data.n_obs, \
        (len(idx_train), len(idx_eval), len(idx_test), test_estim.data.n_obs)
    # 2) Assert that index assignments are exclusive to each split:
    assert len(set(idx_train).intersection(set(idx_eval))) == 0
    assert len(set(idx_train).intersection(set(idx_test))) == 0
    assert len(set(idx_test).intersection(set(idx_eval))) == 0
    # 3) Check partition of index vectors over store data sets matches test split scenario:
    if isinstance(test_estim.estimator.data, DistributedStore):
        # Prepare data set-wise index vectors that are numbered in the same way as global split index vectors.
        # See also EstimatorKeras.train and DistributedStore.subset_cells_idx_global
        idx_raw = test_estim.estimator.data.indices_global
        if isinstance(test_split, float):
            # Make sure that indices from each split are in each data set:
            for z in [idx_train, idx_eval, idx_test]:
                assert np.all([  # in each data set
                    np.any([y in z for y in x])  # at least one match of data set to split index set
                    for x in idx_raw
                ])
        else:
            # Make sure that indices from (train, val) and test split are exclusive:
            datasets_train = np.where([  # in each data set
                np.any([y in idx_train for y in x])  # at least one match of data set to split index set
                for x in idx_raw
            ])[0]
            datasets_eval = np.where([  # in each data set
                np.any([y in idx_eval for y in x])  # at least one match of data set to split index set
                for x in idx_raw
            ])[0]
            datasets_test = np.where([  # in each data set
                np.any([y in idx_test for y in x])  # at least one match of data set to split index set
                for x in idx_raw
            ])[0]
            assert datasets_train == datasets_eval, (datasets_train, datasets_eval)
            assert len(set(datasets_train).intersection(set(datasets_test))) == 0, (datasets_train, datasets_test)
    # 4) Assert that observations mapped to indices are actually unique based on expression vectors:
    # Build numpy arrays of expression input data sets from tensorflow data sets directly from estimator.
    # These data sets are the most processed transformation of the data and stand directly in concat with the model.
    t0 = time.time()
    ds_train = test_estim.estimator._get_dataset(idx=idx_train, batch_size=128, mode='eval', shuffle_buffer_size=1,
                                                 retrieval_batch_size=128)
    print(f"time for building training data set: {time.time() - t0}s")
    t0 = time.time()
    ds_eval = test_estim.estimator._get_dataset(idx=idx_eval, batch_size=128, mode='eval', shuffle_buffer_size=1,
                                                retrieval_batch_size=128)
    print(f"time for building validation data set: {time.time() - t0}s")
    t0 = time.time()
    ds_test = test_estim.estimator._get_dataset(idx=idx_test, batch_size=128, mode='eval', shuffle_buffer_size=1,
                                                retrieval_batch_size=128)
    print(f"time for building test data set: {time.time() - t0}s")
    x_train = []
    x_eval = []
    x_test = []
    t0 = time.time()
    for x, y in ds_train.as_numpy_iterator():
        x_train.append(x[0])
    x_train = np.concatenate(x_train, axis=0)
    print(f"time for iterating over training data set: {time.time() - t0}s")
    t0 = time.time()
    for x, y in ds_eval.as_numpy_iterator():
        x_eval.append(x[0])
    x_eval = np.concatenate(x_eval, axis=0)
    print(f"time for iterating over validation data set: {time.time() - t0}s")
    t0 = time.time()
    for x, y in ds_test.as_numpy_iterator():
        x_test.append(x[0])
    x_test = np.concatenate(x_test, axis=0)
    print(f"time for iterating over test data set: {time.time() - t0}s")
    # Validate size of recovered numpy data sets:
    print(f"shapes received {(x_train.shape[0], x_eval.shape[0], x_test.shape[0])}")
    print(f"shapes expected {(len(idx_train), len(idx_eval), len(idx_test))}")
    assert x_train.shape[0] == len(idx_train)
    assert x_eval.shape[0] == len(idx_eval)
    assert x_test.shape[0] == len(idx_test)
    # Assert that observations are unique within partition:
    assert np.all([
        np.sum([np.all(x_train[i] == x_train[j]) for j in range(x_train.shape[0])]) == 1
        for i in range(x_train.shape[0])
    ])
    assert np.all([
        np.sum([np.all(x_eval[i] == x_eval[j]) for j in range(x_eval.shape[0])]) == 1
        for i in range(x_eval.shape[0])
    ])
    assert np.all([
        np.sum([np.all(x_test[i] == x_test[j]) for j in range(x_test.shape[0])]) == 1
        for i in range(x_test.shape[0])
    ])
    # Assert that observations are not replicated across partitions:
    assert not np.any([
        np.any([np.all(x_train[i] == x_eval[j]) for j in range(x_eval.shape[0])])
        for i in range(x_train.shape[0])
    ])
    assert not np.any([
        np.any([np.all(x_train[i] == x_test[j]) for j in range(x_test.shape[0])])
        for i in range(x_train.shape[0])
    ])
    assert not np.any([
        np.any([np.all(x_test[i] == x_eval[j]) for j in range(x_eval.shape[0])])
        for i in range(x_test.shape[0])
    ])
