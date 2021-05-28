import abc
import anndata
import numpy as np
import os
import pandas as pd
import pytest
import time
from typing import Union

from sfaira.data import load_store, DistributedStoreBase
from sfaira.estimators import EstimatorKeras, EstimatorKerasCelltype, EstimatorKerasEmbedding
from sfaira.versions.genomes import CustomFeatureContainer
from sfaira.versions.metadata import OntologyOboCustom
from sfaira.versions.topologies import TopologyContainer
from sfaira.unit_tests.utils import cached_store_writing, simulate_anndata

dir_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data")
dir_meta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data", "meta")
cache_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
                         "cache", "genomes")

ASSEMBLY = {
    "mouse": "Mus_musculus.GRCm38.102",
    "human": "Homo_sapiens.GRCh38.102",
}
GENES = {
    "mouse": ["ENSMUSG00000000003", "ENSMUSG00000000028"],
    "human": ["ENSG00000000003", "ENSG00000000005"],
}
TARGETS = ["T cell", "CD4-positive helper T cell", "stromal cell", "UNKNOWN"]
TARGET_UNIVERSE = ["CD4-positive helper T cell", "stromal cell"]
ASSAYS = ["10x sequencing", "Smart-seq2"]


TOPOLOGY_EMBEDDING_MODEL = {
    "model_type": None,
    "input": {
        "genome": None,
        "genes": None,
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
        "genome": None,
        "genes": None,
    },
    "output": {
        "cl": "v2021-02-01",
        "targets": TARGET_UNIVERSE
    },
    "hyper_parameters": {
        "l1_coef": 0.,
        "l2_coef": 0.,
    }
}


class HelperEstimatorBase:

    data: Union[anndata.AnnData, DistributedStoreBase]
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

    def load_store(self, organism, organ):
        store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=ASSEMBLY[organism],
                                          organism=organism, organ=organ)
        store = load_store(cache_path=store_path)
        store.subset(attr_key="organism", values=organism)
        store.subset(attr_key="organ", values=organ)
        self.data = store

    @abc.abstractmethod
    def init_topology(self, model_type: str, feature_space: str, organism: str):
        pass

    @abc.abstractmethod
    def init_estimator(self, test_split):
        """
        Initialise target estimator as .estimator attribute.
        """
        pass

    def estimator_train(self, test_split, randomized_batch_access):
        self.estimator.train(
            optimizer="adam",
            lr=0.005,
            epochs=2,
            batch_size=4,
            validation_split=0.5,
            test_split=test_split,
            validation_batch_size=4,
            max_validation_steps=1,
            shuffle_buffer_size=None if randomized_batch_access else 10,
            cache_full=False,
            randomized_batch_access=randomized_batch_access,
        )

    @abc.abstractmethod
    def basic_estimator_test(self, test_split):
        pass

    def load_estimator(self, model_type, data_type, feature_space, test_split, organism="mouse", organ="lung"):
        self.init_topology(model_type=model_type, feature_space=feature_space, organism=organism)
        np.random.seed(1)
        if data_type == "adata":
            self.load_adata()
        else:
            self.load_store(organism=organism, organ=organ)
        self.init_estimator(test_split=test_split)

    def fatal_estimator_test(self, model_type, data_type, test_split=0.1, feature_space="small"):
        self.load_estimator(model_type=model_type, data_type=data_type, feature_space=feature_space,
                            test_split=test_split)
        self.estimator_train(test_split=test_split, randomized_batch_access=False)
        self.basic_estimator_test(test_split=test_split)


class HelperEstimatorKerasEmbedding(HelperEstimatorBase):

    estimator: EstimatorKerasEmbedding
    model_type: str
    tc: TopologyContainer

    def init_topology(self, model_type: str, feature_space: str, organism: str):
        topology = TOPOLOGY_EMBEDDING_MODEL.copy()
        if feature_space == "full":
            # Read 500 genes (not full protein coding) to compromise between being able to distinguish observations
            # and reducing run time of unit tests.
            tab = pd.read_csv(os.path.join(cache_dir, ASSEMBLY[organism] + ".csv"))
            genes_full = tab.loc[tab["gene_biotype"].values == "protein_coding", "gene_id"].values[:500].tolist()
            topology["input"]["genes"] = ["ensg", genes_full]
        else:
            topology["input"]["genes"] = ["ensg", GENES[organism]]
        topology["input"]["genome"] = ASSEMBLY[organism]
        topology["model_type"] = model_type
        if model_type == "linear":
            topology["hyper_parameters"]["latent_dim"] = 2
        else:
            topology["hyper_parameters"]["latent_dim"] = (2, 2, 2)
        self.model_type = model_type
        self.tc = TopologyContainer(topology=topology, topology_id="0.1")

    def init_estimator(self, test_split):
        self.estimator = EstimatorKerasEmbedding(
            data=self.data,
            model_dir=None,
            model_id="testid",
            model_topology=self.tc
        )
        self.estimator.init_model()
        self.estimator.split_train_val_test(test_split=test_split, val_split=0.1)

    def basic_estimator_test(self, test_split=0.1):
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

    def init_topology(self, model_type: str, feature_space: str, organism: str):
        topology = TOPOLOGY_CELLTYPE_MODEL.copy()
        topology["model_type"] = model_type
        topology["input"]["genome"] = ASSEMBLY[organism]
        topology["input"]["genes"] = ["ensg", GENES[organism]]
        if model_type == "mlp":
            topology["hyper_parameters"]["units"] = (2,)
        self.model_type = model_type
        self.tc = TopologyContainer(topology=topology, topology_id="0.0.1")

    def init_estimator(self, test_split):
        tc = self.tc
        if isinstance(self.data, DistributedStoreBase):
            # Reset leaves below:
            tc.topology["output"]["targets"] = None
        self.estimator = EstimatorKerasCelltype(
            data=self.data,
            model_dir=None,
            model_id="testid",
            model_topology=tc
        )
        if isinstance(self.data, DistributedStoreBase):
            leaves = self.estimator.celltype_universe.onto_cl.get_effective_leaves(
                x=[x for x in self.data.obs[self.data._adata_ids_sfaira.cellontology_class].values
                   if x != self.data._adata_ids_sfaira.unknown_celltype_identifier]
            )
            self.nleaves = len(leaves)
            self.estimator.celltype_universe.onto_cl.leaves = leaves
        else:
            self.nleaves = None
        self.estimator.init_model()
        self.estimator.split_train_val_test(test_split=test_split, val_split=0.1)

    def basic_estimator_test(self, test_split=0.1):
        _ = self.estimator.evaluate()
        prediction_output = self.estimator.predict()
        if isinstance(self.estimator.data, anndata.AnnData):
            assert prediction_output.shape[1] == len(TARGET_UNIVERSE), prediction_output.shape
        else:
            assert prediction_output.shape[1] == self.nleaves, prediction_output.shape
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


class HelperEstimatorKerasCelltypeCustomObo(HelperEstimatorKerasCelltype):

    def init_obo_custom(self) -> OntologyOboCustom:
        return OntologyOboCustom(obo=os.path.join(os.path.dirname(__file__), "custom.obo"))

    def init_genome_custom(self, n_features) -> CustomFeatureContainer:
        return CustomFeatureContainer(genome_tab=pd.DataFrame({
            "gene_name": ["dim_" + str(i) for i in range(n_features)],
            "gene_id": ["dim_" + str(i) for i in range(n_features)],
            "gene_biotype": ["embedding" for _ in range(n_features)],
        }))

    def init_topology_custom(self, model_type: str, n_features):
        topology = TOPOLOGY_CELLTYPE_MODEL.copy()
        topology["model_type"] = model_type
        topology["input"]["genome"] = "custom"
        topology["input"]["genes"] = ["biotype", "embedding"]
        topology["output"]["cl"] = "custom"
        topology["output"]["targets"] = ["MYONTO:02", "MYONTO:03"]
        if model_type == "mlp":
            topology["hyper_parameters"]["units"] = (2,)
        self.model_type = model_type
        gc = self.init_genome_custom(n_features=n_features)
        self.tc = TopologyContainer(topology=topology, topology_id="0.0.1", custom_genome_constainer=gc)

    def fatal_estimator_test_custom(self):
        self.init_topology_custom(model_type="mlp", n_features=50)
        obo = self.init_obo_custom()
        np.random.seed(1)
        self.data = simulate_anndata(n_obs=100, genes=self.tc.gc.ensembl,
                                     targets=["MYONTO:01", "MYONTO:02", "MYONTO:03"], obo=obo)
        self.estimator = EstimatorKerasCelltype(
            data=self.data,
            model_dir=None,
            model_id="testid",
            model_topology=self.tc,
            celltype_ontology=obo,
        )
        self.estimator.init_model()
        self.estimator_train(test_split=0.1, randomized_batch_access=False)
        self.basic_estimator_test(test_split=0.1)

# Test embedding models:


@pytest.mark.parametrize("data_type", ["adata", "store"])
def test_for_fatal_linear(data_type):
    test_estim = HelperEstimatorKerasEmbedding()
    test_estim.fatal_estimator_test(model_type="linear", data_type=data_type)


@pytest.mark.parametrize("data_type", ["adata"])
def test_for_fatal_ae(data_type):
    test_estim = HelperEstimatorKerasEmbedding()
    test_estim.fatal_estimator_test(model_type="ae", data_type=data_type)


@pytest.mark.parametrize("data_type", ["adata"])
def test_for_fatal_vae(data_type):
    test_estim = HelperEstimatorKerasEmbedding()
    test_estim.fatal_estimator_test(model_type="vae", data_type=data_type)


# Test cell type predictor models:


@pytest.mark.parametrize("data_type", ["adata", "store"])
def test_for_fatal_mlp(data_type):
    test_estim = HelperEstimatorKerasCelltype()
    test_estim.fatal_estimator_test(model_type="mlp", data_type=data_type)


@pytest.mark.parametrize("data_type", ["adata"])
def test_for_fatal_marker(data_type):
    test_estim = HelperEstimatorKerasCelltype()
    test_estim.fatal_estimator_test(model_type="marker", data_type=data_type)


def test_for_fatal_mlp_custom():
    test_estim = HelperEstimatorKerasCelltypeCustomObo()
    test_estim.fatal_estimator_test_custom()

# Test index sets


@pytest.mark.parametrize("organism", ["human"])
@pytest.mark.parametrize("organ", ["lung"])
@pytest.mark.parametrize("batch_size", [1024, 2048, 4096])
@pytest.mark.parametrize("randomized_batch_access", [False, True])
def test_dataset_size(organism: str, organ: str, batch_size: int, randomized_batch_access: bool):
    """
    Test that tf data set from estimator has same size as generator invoked directly from store based on number of
    observations in emitted batches.

    Tests for batch sizes smaller, equal to and larger than retrieval batch size and with and without randomized
    batch access.
    """
    test_estim = HelperEstimatorKerasEmbedding()
    retrieval_batch_size = 2048
    # Need full feature space here because observations are not necessarily different in small model testing feature
    # space with only two genes:
    test_estim.load_estimator(model_type="linear", data_type="store", feature_space="reduced", test_split=0.2,
                              organism=organism, organ=organ)
    idx_train = test_estim.estimator.idx_train
    shuffle_buffer_size = None if randomized_batch_access else 2
    ds_train = test_estim.estimator._get_dataset(idx=idx_train, batch_size=batch_size, mode='eval',
                                                 shuffle_buffer_size=shuffle_buffer_size,
                                                 retrieval_batch_size=retrieval_batch_size,
                                                 randomized_batch_access=randomized_batch_access)
    x_train_shape = 0
    for x, _ in ds_train.as_numpy_iterator():
        x_train_shape += x[0].shape[0]
    # Define raw store generator on train data to compare and check that it has the same size as tf generator exposed
    # by estimator:
    g_train = test_estim.estimator.data.generator(idx=idx_train, batch_size=retrieval_batch_size,
                                                  randomized_batch_access=randomized_batch_access)
    x_train2_shape = 0
    for x, _ in g_train():
        x_train2_shape += x.shape[0]
    assert x_train_shape == x_train2_shape
    assert x_train_shape == len(idx_train)


@pytest.mark.parametrize("organism", ["mouse"])
@pytest.mark.parametrize("organ", ["lung"])
@pytest.mark.parametrize("data_type", ["adata", "store"])
@pytest.mark.parametrize("randomized_batch_access", [False, True])
@pytest.mark.parametrize("test_split", [0.3, {"assay_sc": "10x sequencing"}])
def test_split_index_sets(organism: str, organ: str, data_type: str, randomized_batch_access: bool, test_split):
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
    test_estim.load_estimator(model_type="linear", data_type=data_type, test_split=test_split, feature_space="full",
                              organism=organism, organ=organ)
    idx_train = test_estim.estimator.idx_train
    idx_eval = test_estim.estimator.idx_eval
    idx_test = test_estim.estimator.idx_test
    print(idx_train)
    print(idx_eval)
    print(idx_test)
    # 1) Assert that index assignment sets sum up to full data set:
    # Make sure that there are no repeated indices in each set.
    assert len(idx_train) == len(np.unique(idx_train))
    assert len(idx_eval) == len(np.unique(idx_eval))
    assert len(idx_test) == len(np.unique(idx_test))
    assert len(idx_train) + len(idx_eval) + len(idx_test) == test_estim.data.n_obs, \
        (len(idx_train), len(idx_eval), len(idx_test), test_estim.data.n_obs)
    if isinstance(test_estim.data, DistributedStoreBase):
        assert np.sum([v.shape[0] for v in test_estim.data.adata_by_key.values()]) == test_estim.data.n_obs
    # 2) Assert that index assignments are exclusive to each split:
    assert len(set(idx_train).intersection(set(idx_eval))) == 0
    assert len(set(idx_train).intersection(set(idx_test))) == 0
    assert len(set(idx_test).intersection(set(idx_eval))) == 0
    # 3) Check partition of index vectors over store data sets matches test split scenario:
    if isinstance(test_estim.estimator.data, DistributedStoreBase):
        # Prepare data set-wise index vectors that are numbered in the same way as global split index vectors.
        # See also EstimatorKeras.train and DistributedStoreBase.subset_cells_idx_global
        idx_raw = test_estim.estimator.data.indices_global.values()
        if isinstance(test_split, float):
            # Make sure that indices from each split are in each data set:
            for i, z in enumerate([idx_train, idx_eval, idx_test]):
                matches = [  # in each data set
                    np.any([y in z for y in x])  # at least one match of data set to split index set
                    for x in idx_raw
                ]
                assert np.all(matches), (i, matches)
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
            assert np.all(datasets_train == datasets_eval), (datasets_train, datasets_eval, datasets_test)
            assert len(set(datasets_train).intersection(set(datasets_test))) == 0, \
                (datasets_train, datasets_eval, datasets_test)
    # 4) Assert that observations mapped to indices are actually unique based on expression vectors:
    # Build numpy arrays of expression input data sets from tensorflow data sets directly from estimator.
    # These data sets are the most processed transformation of the data and stand directly in concat with the model.
    shuffle_buffer_size = None if randomized_batch_access else 2
    t0 = time.time()
    ds_train = test_estim.estimator._get_dataset(idx=idx_train, batch_size=1024, mode='eval',
                                                 shuffle_buffer_size=shuffle_buffer_size,
                                                 retrieval_batch_size=2048,
                                                 randomized_batch_access=randomized_batch_access)
    print(f"time for building training data set: {time.time() - t0}s")
    t0 = time.time()
    ds_eval = test_estim.estimator._get_dataset(idx=idx_eval, batch_size=1024, mode='eval',
                                                shuffle_buffer_size=shuffle_buffer_size,
                                                retrieval_batch_size=2048,
                                                randomized_batch_access=randomized_batch_access)
    print(f"time for building validation data set: {time.time() - t0}s")
    t0 = time.time()
    ds_test = test_estim.estimator._get_dataset(idx=idx_test, batch_size=1024, mode='eval',
                                                shuffle_buffer_size=shuffle_buffer_size,
                                                retrieval_batch_size=2048,
                                                randomized_batch_access=randomized_batch_access)
    print(f"time for building test data set: {time.time() - t0}s")
    # Create two copies of test data set to make sure that re-instantiation of a subset does not cause issues.
    ds_test2 = test_estim.estimator._get_dataset(idx=idx_test, batch_size=1024, mode='eval',
                                                 shuffle_buffer_size=shuffle_buffer_size,
                                                 retrieval_batch_size=2048,
                                                 randomized_batch_access=randomized_batch_access)
    print(f"time for building test data set: {time.time() - t0}s")
    x_train = []
    x_eval = []
    x_test = []
    x_test2_shape = 0
    t0 = time.time()
    for x, _ in ds_train.as_numpy_iterator():
        x_train.append(x[0])
    x_train = np.concatenate(x_train, axis=0)
    print(f"time for iterating over training data set: {time.time() - t0}s")
    t0 = time.time()
    for x, _ in ds_eval.as_numpy_iterator():
        x_eval.append(x[0])
    x_eval = np.concatenate(x_eval, axis=0)
    print(f"time for iterating over validation data set: {time.time() - t0}s")
    t0 = time.time()
    for x, _ in ds_test.as_numpy_iterator():
        x_test.append(x[0])
    x_test = np.concatenate(x_test, axis=0)
    print(f"time for iterating over test data set: {time.time() - t0}s")
    # Assert that duplicate of test data has the same shape:
    for x, _ in ds_test2:
        x_test2_shape += x[0].shape[0]
    assert x_test2_shape == x_test.shape[0]
    # Validate size of recovered numpy data sets:
    print(test_estim.data.n_obs)
    print(f"shapes expected {(len(idx_train), len(idx_eval), len(idx_test))}")
    print(f"shapes received {(x_train.shape[0], x_eval.shape[0], x_test.shape[0])}")
    assert x_train.shape[0] + x_eval.shape[0] + x_test.shape[0] == test_estim.data.n_obs
    assert len(idx_train) + len(idx_eval) + len(idx_test) == test_estim.data.n_obs
    assert x_train.shape[0] == len(idx_train)
    assert x_eval.shape[0] == len(idx_eval)
    assert x_test.shape[0] == len(idx_test)
    # Assert that observations are unique within partition:
    assert np.all([
        np.sum(np.abs(x_train[[i], :] - x_train).sum(axis=1) == 0) == 1
        for i in range(x_train.shape[0])
    ])
    assert np.all([
        np.sum(np.abs(x_eval[[i], :] - x_eval).sum(axis=1) == 0) == 1
        for i in range(x_eval.shape[0])
    ])
    assert np.all([
        np.sum(np.abs(x_test[[i], :] - x_test).sum(axis=1) == 0) == 1
        for i in range(x_test.shape[0])
    ])
    # Assert that observations are not replicated across partitions:
    assert not np.any([
        np.any(np.abs(x_train[[i], :] - x_eval).sum(axis=1) == 0)
        for i in range(x_train.shape[0])
    ])
    assert not np.any([
        np.any(np.abs(x_train[[i], :] - x_test).sum(axis=1) == 0)
        for i in range(x_train.shape[0])
    ])
    assert not np.any([
        np.any(np.abs(x_eval[[i], :] - x_test).sum(axis=1) == 0)
        for i in range(x_eval.shape[0])
    ])
