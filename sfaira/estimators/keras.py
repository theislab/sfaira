import abc
import anndata
import hashlib
import numpy as np
import scipy.sparse
try:
    import tensorflow as tf
except ImportError:
    tf = None
from typing import Union
import os
import warnings
from tqdm import tqdm

from sfaira.consts import AdataIdsSfaira, OCS, AdataIds
from sfaira.data import DistributedStoreBase
from sfaira.models import BasicModelKeras
from sfaira.versions.metadata import CelltypeUniverse, OntologyCl, OntologyObo
from sfaira.versions.topologies import TopologyContainer
from .losses import LossLoglikelihoodNb, LossLoglikelihoodGaussian, LossCrossentropyAgg, KLLoss
from .metrics import custom_mse, custom_negll_nb, custom_negll_gaussian, \
    CustomAccAgg, CustomF1Classwise, CustomFprClasswise, CustomTprClasswise, custom_cce_agg


def prepare_sf(x):
    """
    Uses a minimal size factor of 1e-3 for total counts / 1e4
    """
    if len(x.shape) == 2:
        sf = np.asarray(x.sum(axis=1)).flatten()
    elif len(x.shape) == 1:
        sf = np.asarray(x.sum()).flatten()
    else:
        raise ValueError("x.shape > 2")
    sf = np.log(np.maximum(sf / 1e4, 1e-3))
    return sf


class EstimatorKeras:
    """
    Estimator base class for keras models.
    """
    data: Union[anndata.AnnData, DistributedStoreBase]
    model: Union[BasicModelKeras, None]
    topology_container: Union[TopologyContainer, None]
    model_id: Union[str, None]
    weights: Union[np.ndarray, None]
    model_dir: Union[str, None]
    history: Union[dict, None]
    train_hyperparam: Union[dict, None]
    idx_train: Union[np.ndarray, None]
    idx_eval: Union[np.ndarray, None]
    idx_test: Union[np.ndarray, None]
    adata_ids: AdataIds

    def __init__(
            self,
            data: Union[anndata.AnnData, np.ndarray, DistributedStoreBase],
            model_dir: Union[str, None],
            model_class: str,
            model_id: Union[str, None],
            model_topology: TopologyContainer,
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', ''),
            adata_ids: AdataIds = AdataIdsSfaira()
    ):
        self.data = data
        self.model = None
        self.model_dir = model_dir
        self.model_id = model_id
        self.model_class = model_class
        self.topology_container = model_topology
        # Prepare store with genome container sub-setting:
        if isinstance(self.data, DistributedStoreBase):
            self.data.genome_container = self.topology_container.gc

        self.history = None
        self.train_hyperparam = None
        self.idx_train = None
        self.idx_eval = None
        self.idx_test = None
        self.md5 = weights_md5
        self.cache_path = cache_path
        self._adata_ids = adata_ids

    @property
    def model_type(self):
        return self.topology_container.model_type

    @property
    def organism(self):
        return self.topology_container.organism

    def load_pretrained_weights(self):
        """
        Loads model weights from local directory or zenodo.
        """
        if self.model_dir.startswith('http'):
            # Remote repo
            if not os.path.exists(self.cache_path):
                os.makedirs(self.cache_path)

            import urllib.request
            from urllib.parse import urljoin
            from urllib.error import HTTPError
            try:
                urllib.request.urlretrieve(self.model_dir,
                                           os.path.join(self.cache_path, os.path.basename(self.model_dir))
                                           )
                fn = os.path.join(self.cache_path, os.path.basename(self.model_dir))
            except HTTPError:
                try:
                    urllib.request.urlretrieve(urljoin(self.model_dir, f'{self.model_id}_weights.h5'),
                                               os.path.join(self.cache_path, f'{self.model_id}_weights.h5')
                                               )
                    fn = os.path.join(self.cache_path, f"{self.model_id}_weights.h5")
                except HTTPError:
                    try:
                        urllib.request.urlretrieve(
                            urljoin(self.model_dir, f'{self.model_id}_weights.data-00000-of-00001'),
                            os.path.join(self.cache_path, f'{self.model_id}_weights.data-00000-of-00001')
                        )
                        fn = os.path.join(self.cache_path, f"{self.model_id}_weights.data-00000-of-00001")
                    except HTTPError:
                        raise FileNotFoundError('cannot find remote weightsfile')
        else:
            # Local repo
            if not self.model_dir:
                raise ValueError('the model_id is set but the path to the model is empty')
            if os.path.isfile(self.model_dir) \
                    and not self.model_dir.endswith(".h5") \
                    and not self.model_dir.endswith(".data-00000-of-00001"):
                raise ValueError('weights files saved in h5 format need to have an h5 file extension')

            if os.path.isfile(self.model_dir):
                fn = self.model_dir
            elif os.path.isfile(os.path.join(self.model_dir, f"{self.model_id}_weights.data-00000-of-00001")):
                fn = os.path.join(self.model_dir, f"{self.model_id}_weights.data-00000-of-00001")
            elif os.path.isfile(os.path.join(self.model_dir, f"{self.model_id}_weights.h5")):
                fn = os.path.join(self.model_dir, f"{self.model_id}_weights.h5")
            else:
                raise ValueError('the weightsfile could not be found')

        self._assert_md5_sum(fn, self.md5)
        if fn.endswith(".data-00000-of-00001"):
            self.model.training_model.load_weights(".".join(fn.split(".")[:-1]))
        else:
            self.model.training_model.load_weights(fn)

    def save_weights_to_cache(self):
        if not os.path.exists(os.path.join(self.cache_path, 'weights')):
            os.makedirs(os.path.join(self.cache_path, 'weights'))
        fn = os.path.join(self.cache_path, 'weights', f"{self.model_id}_weights_cache.h5")
        self.model.training_model.save_weights(fn)

    def load_weights_from_cache(self):
        fn = os.path.join(self.cache_path, 'weights', f"{self.model_id}_weights_cache.h5")
        self.model.training_model.load_weights(fn)

    def init_model(self, clear_weight_cache=True, override_hyperpar=None):
        """
        instantiate the model
        :return:
        """
        if clear_weight_cache:
            if os.path.exists(os.path.join(self.cache_path, 'weights')):
                for file in os.listdir(os.path.join(self.cache_path, 'weights')):
                    file_path = os.path.join(os.path.join(self.cache_path, 'weights'), file)
                    os.remove(file_path)

    @staticmethod
    def _assert_md5_sum(
            fn,
            target_md5
    ):
        with open(fn, 'rb') as f:
            hsh = hashlib.md5(f.read()).hexdigest()
        if not hsh == target_md5:
            raise ValueError("md5 of %s did not match expectation" % fn)

    @abc.abstractmethod
    def _get_dataset(
            self,
            idx: Union[np.ndarray, None],
            batch_size: Union[int, None],
            mode: str,
            shuffle_buffer_size: int,
            cache_full: bool,
            weighted: bool,
            retrieval_batch_size: int,
            randomized_batch_access: bool,
            prefetch: Union[int, None],
    ) -> tf.data.Dataset:
        pass

    def _get_class_dict(
            self,
            obs_key: str
    ):
        y = self.data.obs[obs_key]
        for i, val in enumerate(y):
            if type(val) == list:
                y[i] = " / ".join(val)
        labels = np.unique(y)
        label_dict = {}
        for i, label in enumerate(labels):
            label_dict.update({label: float(i)})
        return label_dict

    def _prepare_data_matrix(self, idx: Union[np.ndarray, None]) -> scipy.sparse.csr_matrix:
        """
        Subsets observations x features matrix in .data to observation indices (idx, the split) and features defined
        by topology.

        :param idx: Observation index split.
        :return: Data matrix
        """
        # Check that AnnData is not backed. If backed, assume that these processing steps were done before.
        if self.data.isbacked:
            raise ValueError("tried running backed AnnData object through standard pipeline")

        else:
            # Convert data matrix to csr matrix
            if isinstance(self.data.X, np.ndarray):
                # Change NaN to zero. This occurs for example in concatenation of anndata instances.
                if np.any(np.isnan(self.data.X)):
                    self.data.X[np.isnan(self.data.X)] = 0
                x = scipy.sparse.csr_matrix(self.data.X)
            elif isinstance(self.data.X, scipy.sparse.spmatrix):
                x = self.data.X.tocsr()
            else:
                raise ValueError("data type %s not recognized" % type(self.data.X))

            # Subset cells by provided idx
            if idx is not None:
                x = x[idx, :]

            # If the feature space is already mapped to the right reference, return the data matrix immediately
            if self._adata_ids.mapped_features in self.data.uns_keys() and \
                    self.data.uns[self._adata_ids.mapped_features] == self.topology_container.gc.assembly:
                print(f"found {x.shape[0]} observations")
                return x

            # Compute indices of genes to keep
            data_ids = self.data.var[self._adata_ids.gene_id_ensembl].values.tolist()
            target_ids = self.topology_container.gc.ensembl
            idx_map = np.array([data_ids.index(z) for z in target_ids])
            # Assert that each ID from target IDs appears exactly once in data IDs:
            assert np.all([z in data_ids for z in target_ids]), "not all target feature IDs found in data"
            assert np.all([np.sum(z == np.array(data_ids)) <= 1. for z in target_ids]), \
                "duplicated target feature IDs exist in data"
            # Map feature space.
            x = x[:, idx_map]
            print(f"found {len(idx_map)} intersecting features between {x.shape[1]} features in input data set and"
                  f" {self.topology_container.n_var} features in reference genome")
            print(f"found {x.shape[0]} observations")
            return x

    @abc.abstractmethod
    def _get_loss(self):
        pass

    @abc.abstractmethod
    def _metrics(self):
        pass

    def _compile_models(
            self,
            optimizer: tf.keras.optimizers.Optimizer
    ):
        self.model.training_model.compile(
            optimizer=optimizer,
            loss=self._get_loss(),
            metrics=self._metrics()
        )

    def split_train_val_test(self, val_split: float, test_split: Union[float, dict]):
        # Split training and evaluation data.
        np.random.seed(1)
        all_idx = np.arange(0, self.data.n_obs)  # n_obs is both a property of AnnData and DistributedStoreBase
        if isinstance(test_split, float) or isinstance(test_split, int):
            self.idx_test = np.sort(np.random.choice(
                a=all_idx,
                size=round(self.data.n_obs * test_split),
                replace=False,
            ))
        elif isinstance(test_split, dict):
            in_test = np.ones((self.data.n_obs,), dtype=int) == 1
            for k, v in test_split.items():
                if isinstance(v, bool) or isinstance(v, int) or isinstance(v, list):
                    v = [v]
                if isinstance(self.data, anndata.AnnData):
                    if k not in self.data.obs.columns:
                        raise ValueError(f"Did not find column {k} used to define test set in self.data.")
                    in_test = np.logical_and(in_test, np.array([x in v for x in self.data.obs[k].values]))
                elif isinstance(self.data, DistributedStoreBase):
                    idx = self.data.get_subset_idx_global(attr_key=k, values=v)
                    in_test_k = np.ones((self.data.n_obs,), dtype=int) == 0
                    in_test_k[idx] = True
                    in_test = np.logical_and(in_test, in_test_k)
                else:
                    assert False
            self.idx_test = np.sort(np.where(in_test)[0])
            print(f"Found {len(self.idx_test)} out of {self.data.n_obs} cells that correspond to held out data set")
            print(self.idx_test)
        else:
            raise ValueError("type of test_split %s not recognized" % type(test_split))
        idx_train_eval = np.array([x for x in all_idx if x not in self.idx_test])
        np.random.seed(1)
        self.idx_eval = np.sort(np.random.choice(
            a=idx_train_eval,
            size=round(len(idx_train_eval) * val_split),
            replace=False
        ))
        self.idx_train = np.sort([x for x in idx_train_eval if x not in self.idx_eval])

        # Check that none of the train, test, eval partitions are empty
        if not len(self.idx_test):
            warnings.warn("Test partition is empty!")
        if not len(self.idx_eval):
            raise ValueError("The evaluation dataset is empty.")
        if not len(self.idx_train):
            raise ValueError("The train dataset is empty.")

    def train(
            self,
            optimizer: str,
            lr: float,
            epochs: int = 1000,
            max_steps_per_epoch: Union[int, None] = 20,
            batch_size: int = 128,
            validation_split: float = 0.1,
            test_split: Union[float, dict] = 0.,
            validation_batch_size: int = 256,
            max_validation_steps: Union[int, None] = 10,
            patience: int = 20,
            lr_schedule_min_lr: float = 1e-5,
            lr_schedule_factor: float = 0.2,
            lr_schedule_patience: int = 5,
            shuffle_buffer_size: Union[int, None] = None,
            cache_full: bool = False,
            randomized_batch_access: bool = True,
            retrieval_batch_size: int = 512,
            prefetch: Union[int, None] = 1,
            log_dir: Union[str, None] = None,
            callbacks: Union[list, None] = None,
            weighted: bool = False,
            verbose: int = 2
    ):
        """
        Train model.

        Uses validation loss and maximum number of epochs as termination criteria.

        :param optimizer: str corresponding to tf.keras optimizer to use for fitting.
        :param lr: Learning rate
        :param epochs: refer to tf.keras.models.Model.fit() documentation
        :param max_steps_per_epoch: Maximum steps per epoch.
        :param batch_size: refer to tf.keras.models.Model.fit() documentation
        :param validation_split: refer to tf.keras.models.Model.fit() documentation
            Refers to fraction of training data (1-test_split) to use for validation.
        :param test_split: Fraction of data to set apart for test model before train-validation split.
        :param validation_batch_size: Number of validation data observations to evaluate evaluation metrics on.
        :param max_validation_steps: Maximum number of validation steps to perform.
        :param patience: refer to tf.keras.models.Model.fit() documentation
        :param lr_schedule_min_lr: Minimum learning rate for learning rate reduction schedule.
        :param lr_schedule_factor: Factor to reduce learning rate by within learning rate reduction schedule
            when plateau is reached.
        :param lr_schedule_patience: Patience for learning rate reduction in learning rate reduction schedule.
        :param shuffle_buffer_size: tf.Dataset.shuffle(): buffer_size argument.
        :param cache_full: Whether to use tensorflow caching on full training and validation data.
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
        :param log_dir: Directory to save tensorboard callback to. Disabled if None.
        :param callbacks: Add additional callbacks to the training call
        :param weighted:
        :param verbose:
        :return:
        """
        # Set optimizer
        if optimizer.lower() == "adam":
            optim = tf.keras.optimizers.Adam(learning_rate=lr)
        elif optimizer.lower() == "sgd":
            optim = tf.keras.optimizers.SGD(learning_rate=lr)
        elif optimizer.lower() == "rmsprop":
            optim = tf.keras.optimizers.RMSprop(learning_rate=lr)
        elif optimizer.lower() == "adagrad":
            optim = tf.keras.optimizers.Adagrad(learning_rate=lr)
        else:
            assert False
        # Save training settings to allow model restoring.
        self.train_hyperparam = {
            "epochs": epochs,
            "max_steps_per_epoch": max_steps_per_epoch,
            "optimizer": optimizer,
            "lr": lr,
            "batch_size": batch_size,
            "validation_split": validation_split,
            "validation_batch_size": validation_batch_size,
            "max_validation_steps": max_validation_steps,
            "patience": patience,
            "lr_schedule_min_lr": lr_schedule_min_lr,
            "lr_schedule_factor": lr_schedule_factor,
            "lr_schedule_patience": lr_schedule_patience,
            "log_dir": log_dir,
            "weighted": weighted
        }

        # Set callbacks.
        cbs = [tf.keras.callbacks.TerminateOnNaN()]
        if patience is not None and patience > 0:
            cbs.append(tf.keras.callbacks.EarlyStopping(
                monitor='val_loss',
                patience=patience,
                restore_best_weights=True,
                verbose=verbose
            ))
        if lr_schedule_factor is not None and lr_schedule_factor < 1.:
            cbs.append(tf.keras.callbacks.ReduceLROnPlateau(
                monitor='val_loss',
                factor=lr_schedule_factor,
                patience=lr_schedule_patience,
                min_lr=lr_schedule_min_lr,
                verbose=verbose
            ))
        if log_dir is not None:
            cbs.append(tf.keras.callbacks.TensorBoard(
                log_dir=log_dir,
                histogram_freq=0,
                batch_size=32,
                write_graph=False,
                write_grads=False,
                write_images=False,
                embeddings_freq=0,
                embeddings_layer_names=None,
                embeddings_metadata=None,
                embeddings_data=None,
                update_freq='epoch'
            ))

        if callbacks is not None:
            # callbacks needs to be a list
            cbs += callbacks

        # Check randomisation settings:
        if shuffle_buffer_size is not None and shuffle_buffer_size > 0 and randomized_batch_access:
            raise ValueError("You are using shuffle_buffer_size and randomized_batch_access, this is likely not "
                             "intended.")
        if cache_full and randomized_batch_access:
            raise ValueError("You are using cache_full and randomized_batch_access, this is likely not intended.")
        self.split_train_val_test(val_split=validation_split, test_split=test_split)
        self._compile_models(optimizer=optim)
        shuffle_buffer_size = shuffle_buffer_size if shuffle_buffer_size is not None else 0
        train_dataset = self._get_dataset(
            idx=self.idx_train,
            batch_size=batch_size,
            retrieval_batch_size=retrieval_batch_size,
            mode='train',
            shuffle_buffer_size=min(shuffle_buffer_size, len(self.idx_train)),
            weighted=weighted,
            cache_full=cache_full,
            randomized_batch_access=randomized_batch_access,
            prefetch=prefetch,
        )
        eval_dataset = self._get_dataset(
            idx=self.idx_eval,
            batch_size=validation_batch_size,
            retrieval_batch_size=retrieval_batch_size,
            mode='train_val',
            shuffle_buffer_size=min(shuffle_buffer_size, len(self.idx_eval)),
            weighted=weighted,
            cache_full=cache_full,
            randomized_batch_access=randomized_batch_access,
            prefetch=prefetch,
        )

        steps_per_epoch = min(max(len(self.idx_train) // batch_size, 1), max_steps_per_epoch)
        validation_steps = min(max(len(self.idx_eval) // validation_batch_size, 1), max_validation_steps)

        self.history = self.model.training_model.fit(
            x=train_dataset,
            epochs=epochs,
            steps_per_epoch=steps_per_epoch,
            callbacks=cbs,
            validation_data=eval_dataset,
            validation_steps=validation_steps,
            verbose=verbose
        ).history

    @property
    def using_store(self) -> bool:
        return isinstance(self.data, DistributedStoreBase)

    @property
    def obs_train(self):
        return self.data.obs.iloc[self.idx_train, :]

    @property
    def obs_eval(self):
        return self.data.obs.iloc[self.idx_eval, :]

    @property
    def obs_test(self):
        return self.data.obs.iloc[self.idx_test, :]


class EstimatorKerasEmbedding(EstimatorKeras):
    """
    Estimator class for the embedding model.
    """

    def __init__(
            self,
            data: Union[anndata.AnnData, np.ndarray, DistributedStoreBase],
            model_dir: Union[str, None],
            model_id: Union[str, None],
            model_topology: TopologyContainer,
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', ''),
            adata_ids: AdataIds = AdataIdsSfaira()
    ):
        super(EstimatorKerasEmbedding, self).__init__(
            data=data,
            model_dir=model_dir,
            model_class="embedding",
            model_id=model_id,
            model_topology=model_topology,
            weights_md5=weights_md5,
            cache_path=cache_path,
            adata_ids=adata_ids
        )

    def init_model(
            self,
            clear_weight_cache: bool = True,
            override_hyperpar: Union[dict, None] = None
    ):
        """
        instantiate the model
        :return:
        """
        super().init_model(clear_weight_cache=clear_weight_cache)
        if self.model_type == 'vae':
            from sfaira.models.embedding import ModelVaeVersioned as Model
        elif self.model_type == 'ae':
            from sfaira.models.embedding import ModelAeVersioned as Model
        elif self.model_type == 'linear' or self.model_type == 'nmf':
            from sfaira.models.embedding import ModelLinearVersioned as Model
        elif self.model_type == 'vaeiaf':
            from sfaira.models.embedding import ModelVaeIAFVersioned as Model
        elif self.model_type == 'vaevamp':
            from sfaira.models.embedding import ModelVaeVampVersioned as Model
        else:
            raise ValueError(f'unknown model type {self.model_type} for EstimatorKerasEmbedding')
        self.model = Model(
            topology_container=self.topology_container,
            override_hyperpar=override_hyperpar
        )

    @staticmethod
    def _get_output_dim(n_features, model_type, mode='train'):
        if mode == 'predict':  # Output shape is same for predict mode regardless of model type
            output_types = (tf.float32, tf.float32),
            output_shapes = (n_features, ()),
        elif model_type == "vae":
            output_types = ((tf.float32, tf.float32), (tf.float32, tf.float32))
            output_shapes = ((n_features, ()), (n_features, ()))
        else:
            output_types = ((tf.float32, tf.float32), tf.float32)
            output_shapes = ((n_features, ()), n_features)

        return output_types, output_shapes

    def _get_base_generator(
            self,
            generator_helper,
            idx: Union[np.ndarray, None],
            batch_size: int,
            randomized_batch_access: bool,
    ):
        """
        Yield a basic generator based on which a tf dataset can be built.

        The signature of this generator can be modified through generator_helper.

        :param generator_helper: Python function that should take (x_sample,) as an input:

            - x_sample is a gene expression vector of a cell
        :param idx: Indicies of data set to include in generator.
        :param batch_size: Number of observations read from disk in each batched access.
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
        :return:
        """
        if idx is None:
            idx = np.arange(0, self.data.n_obs)

        # Prepare data reading according to whether anndata is backed or not:
        if self.using_store:
            generator_raw = self.data.generator(
                idx=idx,
                batch_size=batch_size,
                obs_keys=[],
                return_dense=True,
                randomized_batch_access=randomized_batch_access,
            )

            def generator():
                for z in generator_raw():
                    x_sample = z[0]
                    if isinstance(x_sample, scipy.sparse.csr_matrix):
                        x_sample = x_sample.todense()
                    x_sample = np.asarray(x_sample)
                    for i in range(x_sample.shape[0]):
                        yield generator_helper(x_sample=x_sample[i])

            n_features = self.data.n_vars
            n_samples = self.data.n_obs
        else:
            x = self.data.X if self.data.isbacked else self._prepare_data_matrix(idx=idx)
            indices = idx if self.data.isbacked else range(x.shape[0])
            n_obs = len(indices)
            remainder = n_obs % batch_size
            batch_starts_ends = [
                (int(x * batch_size), int(x * batch_size) + batch_size)
                for x in np.arange(0, n_obs // batch_size + int(remainder > 0))
            ]

            def generator():
                is_sparse = isinstance(x[0, :], scipy.sparse.spmatrix)
                for s, e in batch_starts_ends:
                    x_sample = np.asarray(x[indices[s:e], :].todense()) if is_sparse \
                        else x[indices[s:e], :]
                    for i in range(x_sample.shape[0]):
                        yield generator_helper(x_sample=x_sample[i])

            n_features = x.shape[1]
            n_samples = x.shape[0]

        return generator, n_samples, n_features

    def _get_dataset(
            self,
            idx: Union[np.ndarray, None],
            batch_size: Union[int, None],
            mode: str,
            shuffle_buffer_size: int = int(1e7),
            cache_full: bool = False,
            weighted: bool = False,
            retrieval_batch_size: int = 128,
            randomized_batch_access: bool = False,
            prefetch: Union[int, None] = 1,
    ) -> tf.data.Dataset:
        """

        :param idx:
        :param batch_size:
        :param mode:
        :param shuffle_buffer_size:
        :param weighted: Whether to use weights. Not implemented for embedding models yet.
        :param retrieval_batch_size: Number of observations read from disk in each batched access.
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
        :return:
        """
        # Determine model type [ae, vae(iaf, vamp)]
        model_type = "vae" if self.model_type[:3] == "vae" else "ae"

        if mode in ['train', 'train_val', 'eval', 'predict']:
            def generator_helper(x_sample):
                sf_sample = prepare_sf(x=x_sample)[0]
                if mode == 'predict':
                    return (x_sample, sf_sample),
                elif model_type == "vae":
                    return (x_sample, sf_sample), (x_sample, sf_sample)
                else:
                    return (x_sample, sf_sample), x_sample

            generator, n_samples, n_features = self._get_base_generator(
                generator_helper=generator_helper,
                idx=idx,
                batch_size=retrieval_batch_size,
                randomized_batch_access=randomized_batch_access,
            )
            output_types, output_shapes = self._get_output_dim(n_features=n_features, model_type=model_type, mode=mode)
            dataset = tf.data.Dataset.from_generator(
                generator=generator,
                output_types=output_types,
                output_shapes=output_shapes
            )
            if cache_full:
                dataset = dataset.cache()
            # Only shuffle in train modes
            if mode in ['train', 'train_val']:
                dataset = dataset.repeat()
                if shuffle_buffer_size is not None and shuffle_buffer_size > 0:
                    dataset = dataset.shuffle(
                        buffer_size=min(n_samples, shuffle_buffer_size),
                        seed=None,
                        reshuffle_each_iteration=True)
            if prefetch is None:
                prefetch = tf.data.AUTOTUNE
            dataset = dataset.batch(batch_size, drop_remainder=False).prefetch(prefetch)

            return dataset

        elif mode == 'gradient_method':  # TODO depreceate this code
            # Prepare data reading according to whether anndata is backed or not:
            cell_to_class = self._get_class_dict(obs_key=self._adata_ids.cellontology_class)
            if self.using_store:
                n_features = self.data.n_vars
                generator_raw = self.data.generator(
                    idx=idx,
                    batch_size=1,
                    obs_keys=[self._adata_ids.cellontology_class],
                    return_dense=True,
                )

                def generator():
                    for z in generator_raw():
                        x_sample = z[0]
                        if isinstance(x_sample, scipy.sparse.csr_matrix):
                            x_sample = x_sample.todense()
                        x_sample = np.asarray(x_sample).flatten()
                        sf_sample = prepare_sf(x=x_sample)[0]
                        y_sample = z[1][self._adata_ids.cellontology_class].values[0]
                        yield (x_sample, sf_sample), (x_sample, cell_to_class[y_sample])

            elif isinstance(self.data, anndata.AnnData) and self.data.isbacked:
                if idx is None:
                    idx = np.arange(0, self.data.n_obs)
                n_features = self.data.X.shape[1]

                def generator():
                    sparse = isinstance(self.data.X[0, :], scipy.sparse.spmatrix)
                    for i in idx:
                        x_sample = self.data.X[i, :].toarray().flatten() if sparse else self.data.X[i, :].flatten()
                        sf_sample = prepare_sf(x=x_sample)[0]
                        y_sample = self.data.obs[self._adata_ids.cellontology_id][i]
                        yield (x_sample, sf_sample), (x_sample, cell_to_class[y_sample])
            else:
                if idx is None:
                    idx = np.arange(0, self.data.n_obs)
                x = self._prepare_data_matrix(idx=idx)
                sf = prepare_sf(x=x)
                y = self.data.obs[self._adata_ids.cellontology_class].values[idx]
                # for gradients per celltype in compute_gradients_input()
                n_features = x.shape[1]

                def generator():
                    for i in range(x.shape[0]):
                        yield (x[i, :].toarray().flatten(), sf[i]), (x[i, :].toarray().flatten(), cell_to_class[y[i]])

            output_types, output_shapes = self._get_output_dim(n_features, 'vae')
            dataset = tf.data.Dataset.from_generator(
                generator=generator,
                output_types=output_types,
                output_shapes=output_shapes
            )
            dataset = dataset.shuffle(
                buffer_size=shuffle_buffer_size,
                seed=None,
                reshuffle_each_iteration=True
            ).batch(batch_size, drop_remainder=False).prefetch(prefetch)

            return dataset

        else:
            raise ValueError(f'Mode {mode} not recognised. Should be "train", "eval" or" predict"')

    def _get_loss(self):
        if self.topology_container.topology["hyper_parameters"]["output_layer"] in [
            "nb", "nb_const_disp", "nb_shared_disp"
        ]:
            reconstruction_loss_vae = LossLoglikelihoodNb(average=False)  # TODO maybe handly kwargs via functional?
            reconstruction_loss = LossLoglikelihoodNb(average=True)
        elif self.topology_container.topology["hyper_parameters"]["output_layer"] in [
            "gaussian", "gaussian_const_disp", "gaussian_shared_disp"
        ]:
            reconstruction_loss_vae = LossLoglikelihoodGaussian(average=False)
            reconstruction_loss = LossLoglikelihoodGaussian(average=True)
        else:
            raise ValueError(self.topology_container.topology["hyper_parameters"]["output_layer"])

        if self.model_type[:3] == "vae":  # TODO too hacky
            return {
                "neg_ll": reconstruction_loss_vae,
                "kl": KLLoss()
            }
        else:
            return {"neg_ll": reconstruction_loss}

    def _metrics(self):
        if self.topology_container.topology["hyper_parameters"]["output_layer"] in [
            "nb", "nb_const_disp", "nb_shared_disp"
        ]:
            custom_negll = custom_negll_nb
        elif self.topology_container.topology["hyper_parameters"]["output_layer"] in [
            "gaussian", "gaussian_const_disp", "gaussian_shared_disp"
        ]:
            custom_negll = custom_negll_gaussian
        else:
            raise ValueError(self.topology_container.topology["hyper_parameters"]["output_layer"])

        custom_negll.__name__ = "custom_negll"

        return {"neg_ll": [custom_mse, custom_negll]}

    def evaluate_any(self, idx, batch_size: int = 128, max_steps: int = np.inf):
        """
        Evaluate the custom model on any local data.

        Defaults to run on full data if idx is None.

        :param idx: Indices of observations to evaluate on. Evaluates on all observations if None.
        :param batch_size: Batch size for evaluation.
        :param max_steps: Maximum steps before evaluation round is considered complete.
        :return: Dictionary of metric names and values.
        """
        if idx is None or idx.any():  # true if the array is not empty or if the passed value is None
            idx = np.arange(0, self.data.n_obs) if idx is None else idx
            dataset = self._get_dataset(
                idx=idx,
                batch_size=batch_size,
                mode='eval',
                retrieval_batch_size=128,
                shuffle_buffer_size=0,
            )
            steps = min(max(len(idx) // batch_size, 1), max_steps)
            results = self.model.training_model.evaluate(x=dataset, steps=steps)
            return dict(zip(self.model.training_model.metrics_names, results))
        else:
            return {}

    def evaluate(self, batch_size: int = 128, max_steps: int = np.inf):
        """
        Evaluate the custom model on test data.

        Defaults to run on full data if idx_test was not set before, ie. train() has not been called before.

        :param batch_size: Batch size for evaluation.
        :param max_steps: Maximum steps before evaluation round is considered complete.
        :return: Dictionary of metric names and values.
        """
        return self.evaluate_any(idx=self.idx_test, batch_size=batch_size, max_steps=max_steps)

    def predict(self, batch_size: int = 128):
        """
        return the prediction of the model

        :return:
        prediction
        """
        if self.idx_test is None or self.idx_test.any():
            dataset = self._get_dataset(
                idx=self.idx_test,
                batch_size=batch_size,
                mode='predict',
                retrieval_batch_size=128,
                shuffle_buffer_size=0,
            )
            return self.model.predict_reconstructed(x=dataset)
        else:
            return np.array([])

    def predict_embedding(self, batch_size: int = 128):
        """
        return the prediction in the latent space (z_mean for variational models)

        :return:
        latent space
        """
        if self.idx_test is None or self.idx_test.any():
            dataset = self._get_dataset(
                idx=self.idx_test,
                batch_size=batch_size,
                mode='predict',
                retrieval_batch_size=128,
                shuffle_buffer_size=0,
            )
            return self.model.predict_embedding(x=dataset, variational=False)
        else:
            return np.array([])

    def predict_embedding_variational(self, batch_size: int = 128, max_steps: int = np.inf):
        """
        return the prediction of z, z_mean, z_log_var in the variational latent space

        :return:
        sample of latent space, mean of latent space, variance of latent space
        """
        if self.idx_test is None or self.idx_test:
            dataset = self._get_dataset(
                idx=self.idx_test,
                batch_size=batch_size,
                mode='predict',
                retrieval_batch_size=128,
                shuffle_buffer_size=0,
            )
            return self.model.predict_embedding(x=dataset, variational=True)
        else:
            return np.array([])

    def compute_gradients_input(
            self,
            batch_size: int = 128,
            test_data: bool = False,
            abs_gradients: bool = True,
            per_celltype: bool = False
    ):
        if test_data:
            idx = self.idx_test
            if self.idx_test is None:
                num_samples = 10000
                idx = np.random.randint(0, self.data.X.shape[0], num_samples)
            n_obs = len(idx)
        else:
            idx = None
            n_obs = self.data.X.shape[0]

        ds = self._get_dataset(
            idx=idx,
            batch_size=batch_size,
            mode='gradient_method',  # to get a tf.GradientTape compatible data set
        )

        if per_celltype:
            cell_to_id = self._get_class_dict(obs_key=self._adata_ids.cellontology_class)
            cell_names = cell_to_id.keys()
            cell_id = cell_to_id.values()
            id_to_cell = dict([(key, value) for (key, value) in zip(cell_id, cell_names)])
            grads_x = dict([(key, 0) for key in cell_names])
            counts = dict([(key, 0) for key in cell_names])
        else:
            grads_x = 0
        # Loop over sub-selected data set and sum gradients across all selected observations.
        if self.model_type[:3] == "vae":  # TODO: fix bug for vaeiaf model. This function can not be called for vaeiaf model
            model = tf.keras.Model(
                self.model.training_model.input,
                self.model.encoder_model.output[0]
            )
            latent_dim = self.model.encoder_model.output[0].shape[1]
            input_dim = self.model.training_model.input[0].shape[1]
        else:
            model = tf.keras.Model(
                self.model.training_model.input,
                self.model.encoder_model.output
            )
            latent_dim = self.model.encoder_model.output[0].shape[0]
            input_dim = self.model.training_model.input[0].shape[1]

        @tf.function
        def get_gradients(x_batch):
            x, sf = x_batch
            with tf.GradientTape(persistent=True) as tape:
                tape.watch(x)
                model_out = model((x, sf))
            if abs_gradients:
                def f(xx):
                    return abs(xx)
            else:
                def f(xx):
                    return xx
            # marginalize on batch level and then accumulate batches
            # batch_jacobian gives output of size: (batch_size, latent_dim, input_dim)
            return f(tape.batch_jacobian(model_out, x))

        for step, (x_batch, y_batch) in tqdm(enumerate(ds), total=np.ceil(n_obs / batch_size)):
            batch_gradients = get_gradients(x_batch).numpy()
            _, y = y_batch
            if per_celltype:
                for i in np.unique(y):
                    grads_x[id_to_cell[i]] += np.sum(batch_gradients[y == i, :, :], axis=0)
                    counts[id_to_cell[i]] += np.sum(y == i)
            else:
                grads_x += np.sum(batch_gradients, axis=0)
        if per_celltype:
            for cell in cell_names:
                print(f'{cell} with {counts[cell]} observations')
                grads_x[cell] = grads_x[cell] / counts[cell] if counts[cell] > 0 else np.zeros((latent_dim, input_dim))

            return {'gradients': grads_x, 'counts': counts}
        else:
            return grads_x / n_obs


class EstimatorKerasCelltype(EstimatorKeras):
    """
    Estimator class for the cell type model.
    """

    celltype_universe: CelltypeUniverse

    def __init__(
            self,
            data: Union[anndata.AnnData, DistributedStoreBase],
            model_dir: Union[str, None],
            model_id: Union[str, None],
            model_topology: TopologyContainer,
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', ''),
            celltype_ontology: Union[OntologyObo, None] = None,
            max_class_weight: float = 1e3,
            remove_unlabeled_cells: bool = True,
            adata_ids: AdataIds = AdataIdsSfaira()
    ):
        super(EstimatorKerasCelltype, self).__init__(
            data=data,
            model_dir=model_dir,
            model_class="celltype",
            model_id=model_id,
            model_topology=model_topology,
            weights_md5=weights_md5,
            cache_path=cache_path,
            adata_ids=adata_ids
        )
        if remove_unlabeled_cells:
            # Remove cells without type label from store:
            if isinstance(self.data, DistributedStoreBase):
                self.data.subset(attr_key="cellontology_class", excluded_values=[
                    self._adata_ids.unknown_celltype_identifier,
                    self._adata_ids.not_a_cell_celltype_identifier,
                    None,  # TODO: it may be possible to remove this in the future
                    np.nan,  # TODO: it may be possible to remove this in the future
                ])
            elif isinstance(self.data, anndata.AnnData):
                self.data = self.data[np.where([
                    x not in [
                        self._adata_ids.unknown_celltype_identifier,
                        self._adata_ids.not_a_cell_celltype_identifier,
                        None,  # TODO: it may be possible to remove this in the future
                        np.nan,  # TODO: it may be possible to remove this in the future
                    ] for x in self.data.obs[self._adata_ids.cellontology_class].values
                ])[0], :]
            else:
                assert False
        assert "cl" in self.topology_container.output.keys(), self.topology_container.output.keys()
        assert "targets" in self.topology_container.output.keys(), self.topology_container.output.keys()
        self.max_class_weight = max_class_weight
        if celltype_ontology is None:
            celltype_ontology = OntologyCl(branch=self.topology_container.output["cl"])
        self.celltype_universe = CelltypeUniverse(
            cl=celltype_ontology,
            uberon=OCS.organ,
        )
        # Set leaves if they are defined in topology:
        if self.topology_container.output["targets"] is not None:
            self.celltype_universe.onto_cl.leaves = self.topology_container.output["targets"]

    def init_model(
            self,
            clear_weight_cache: bool = True,
            override_hyperpar: Union[dict, None] = None
    ):
        """
        instantiate the model
        :return:
        """
        super().init_model(clear_weight_cache=clear_weight_cache)
        if self.model_type == "marker":
            from sfaira.models.celltype import CellTypeMarkerVersioned as Model
        elif self.model_type == "mlp":
            from sfaira.models.celltype import CellTypeMlpVersioned as Model
        else:
            raise ValueError('unknown topology %s for EstimatorKerasCelltype' % self.model_type)

        self.model = Model(
            celltypes_version=self.celltype_universe,
            topology_container=self.topology_container,
            override_hyperpar=override_hyperpar
        )

    @property
    def ntypes(self):
        return self.celltype_universe.onto_cl.n_leaves

    @property
    def ontology_ids(self):
        return self.celltype_universe.onto_cl.convert_to_id(self.celltype_universe.onto_cl.leaves)

    @property
    def ontology_names(self):
        return self.celltype_universe.onto_cl.convert_to_name(self.celltype_universe.onto_cl.leaves)

    def _one_hot_encoder(self):
        leave_maps = self.celltype_universe.onto_cl.prepare_maps_to_leaves(include_self=True)

        def encoder(x) -> np.ndarray:
            if isinstance(x, str):
                x = [x]
            # Encodes unknowns to empty rows.
            idx = [
                leave_maps[y] if y not in [
                    self._adata_ids.unknown_celltype_identifier,
                    self._adata_ids.not_a_cell_celltype_identifier,
                ] else np.array([])
                for y in x
            ]
            oh = np.zeros((len(x), self.ntypes,), dtype="float32")
            for i, y in enumerate(idx):
                scale = len(y)
                if scale > 0:
                    oh[i, y] = 1. / scale
            return oh

        return encoder

    def _get_celltype_out(
            self,
            idx: Union[np.ndarray, None],
    ):
        """
        Build one hot encoded cell type output tensor and observation-wise weight matrix.

        :param lookup_ontology: list of ontology names to consider.
        :return:
        """
        if idx is None:
            idx = np.arange(0, self.data.n_obs)
        # One whether "unknown" is already included, otherwise add one extra column.
        onehot_encoder = self._one_hot_encoder()
        y = np.concatenate([
            onehot_encoder(z)
            for z in self.data.obs[self._adata_ids.cellontology_id].values[idx].tolist()
        ], axis=0)
        # Distribute aggregated class weight for computation of weights:
        freq = np.mean(y / np.sum(y, axis=1, keepdims=True), axis=0, keepdims=True)
        weights = 1. / np.matmul(y, freq.T)  # observation wise weight matrix
        # Threshold weights:
        weights = np.asarray(
            np.minimum(weights, np.zeros_like(weights) + self.max_class_weight),
            dtype="float32"
        ).flatten()
        return weights, y

    @staticmethod
    def _get_output_dim(n_features, n_labels, mode):
        if mode == 'predict':
            output_types = (tf.float32,)
            output_shapes = (tf.TensorShape([n_features]),)
        else:
            output_types = (tf.float32, tf.float32, tf.float32)
            output_shapes = (
                (tf.TensorShape([n_features])),
                tf.TensorShape([n_labels]),
                tf.TensorShape([])
            )

        return output_types, output_shapes

    def _get_base_generator(
            self,
            idx: Union[np.ndarray, None],
            yield_labels: bool,
            weighted: bool,
            batch_size: int,
            randomized_batch_access: bool,
    ):
        """
        Yield a basic generator based on which a tf dataset can be built.

        The signature of this generator can be modified through generator_helper.

        :param generator_helper: Python function that should take (x_sample, y_sample, w_sample) as an input:

            - x_sample is a gene expression vector of a cell
            - y_sample is a one-hot encoded label vector of a cell
            - w_sample is a weight scalar of a cell
        :param idx: Indicies of data set to include in generator.
        :param yield_labels:
        :param batch_size: Number of observations read from disk in each batched access.
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
        :return:
        """
        if idx is None:
            idx = np.arange(0, self.data.n_obs)

        # Prepare data reading according to whether anndata is backed or not:
        if self.using_store:
            if weighted:
                raise ValueError("using weights with store is not supported yet")
            generator_raw = self.data.generator(
                idx=idx,
                batch_size=batch_size,
                obs_keys=[self._adata_ids.cellontology_id],
                return_dense=True,
                randomized_batch_access=randomized_batch_access,
            )
            if yield_labels:
                onehot_encoder = self._one_hot_encoder()

            def generator():
                for z in generator_raw():
                    x_sample = z[0]
                    if isinstance(x_sample, scipy.sparse.csr_matrix):
                        x_sample = x_sample.todense()
                    x_sample = np.asarray(x_sample)
                    if yield_labels:
                        y_sample = onehot_encoder(z[1][self._adata_ids.cellontology_id].values)
                        for i in range(x_sample.shape[0]):
                            if y_sample[i].sum() > 0:
                                yield x_sample[i], y_sample[i], 1.
                    else:
                        for i in range(x_sample.shape[0]):
                            yield x_sample[i],
            n_features = self.data.n_vars
            n_samples = self.data.n_obs
        else:
            if yield_labels:
                weights, y = self._get_celltype_out(idx=idx)
                if not weighted:
                    weights = np.ones_like(weights)
            x = self.data.X if self.data.isbacked else self._prepare_data_matrix(idx=idx)
            is_sparse = isinstance(x, scipy.sparse.spmatrix)
            indices = idx if self.data.isbacked else range(x.shape[0])
            n_obs = len(indices)
            remainder = n_obs % batch_size
            batch_starts_ends = [
                (int(x * batch_size), int(x * batch_size) + batch_size)
                for x in np.arange(0, n_obs // batch_size + int(remainder > 0))
            ]

            def generator():
                for s, e in batch_starts_ends:
                    x_sample = np.asarray(x[indices[s:e], :].todense()) if is_sparse else x[indices[s:e], :]
                    if yield_labels:
                        y_sample = y[indices[s:e], :]
                        w_sample = weights[indices[s:e]]
                        for i in range(x_sample.shape[0]):
                            if y_sample[i].sum() > 0:
                                yield x_sample[i], y_sample[i], w_sample[i]
                    else:
                        for i in range(x_sample.shape[0]):
                            yield x_sample[i],

            n_features = x.shape[1]
            n_samples = x.shape[0]

        n_labels = self.celltype_universe.onto_cl.n_leaves
        return generator, n_samples, n_features, n_labels

    def _get_dataset(
            self,
            idx: Union[np.ndarray, None],
            batch_size: Union[int, None],
            mode: str,
            shuffle_buffer_size: int = int(1e7),
            cache_full: bool = False,
            weighted: bool = False,
            retrieval_batch_size: int = 128,
            randomized_batch_access: bool = False,
            prefetch: Union[int, None] = 1,
    ) -> tf.data.Dataset:
        """

        :param idx:
        :param batch_size:
        :param mode:
        :param shuffle_buffer_size:
        :param weighted: Whether to use weights.
        :param retrieval_batch_size: Number of observations read from disk in each batched access.
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
        :return:
        """
        generator, n_samples, n_features, n_labels = self._get_base_generator(
            idx=idx,
            yield_labels=mode in ['train', 'train_val', 'eval'],
            weighted=weighted,
            batch_size=retrieval_batch_size,
            randomized_batch_access=randomized_batch_access,
        )
        output_types, output_shapes = self._get_output_dim(n_features=n_features, n_labels=n_labels, mode=mode)
        dataset = tf.data.Dataset.from_generator(
            generator=generator,
            output_types=output_types,
            output_shapes=output_shapes
        )
        if cache_full:
            dataset = dataset.cache()
        if mode == 'train' or mode == 'train_val':
            dataset = dataset.repeat()
            if shuffle_buffer_size is not None and shuffle_buffer_size > 0:
                dataset = dataset.shuffle(
                    buffer_size=min(n_samples, shuffle_buffer_size),
                    seed=None,
                    reshuffle_each_iteration=True)
        if prefetch is None:
            prefetch = tf.data.AUTOTUNE
        dataset = dataset.batch(batch_size, drop_remainder=False).prefetch(prefetch)

        return dataset

    def _get_loss(self):
        return LossCrossentropyAgg()

    def _metrics(self):
        return [
            "accuracy",
            custom_cce_agg,
            CustomAccAgg(),
            CustomF1Classwise(k=self.ntypes),
            CustomFprClasswise(k=self.ntypes),
            CustomTprClasswise(k=self.ntypes)
        ]

    def predict(self, batch_size: int = 128, max_steps: int = np.inf):
        """
        Return the prediction of the model

        :param batch_size: Batch size for evaluation.
        :param max_steps: Maximum steps before evaluation round is considered complete.
        :return: Prediction tensor.
        """
        idx = self.idx_test
        if idx is None or idx.any():
            dataset = self._get_dataset(
                idx=idx,
                batch_size=batch_size,
                mode='predict',
                retrieval_batch_size=128,
                shuffle_buffer_size=0,
            )
            return self.model.training_model.predict(x=dataset)
        else:
            return np.array([])

    def ytrue(self, batch_size: int = 128, max_steps: int = np.inf):
        """
        Return the true labels of the test set.

        :return: true labels
        """
        if self.idx_test is None or self.idx_test.any():
            dataset = self._get_dataset(
                idx=self.idx_test,
                batch_size=batch_size,
                mode='eval',
                shuffle_buffer_size=0,
            )
            y_true = []
            for _, y, _ in dataset.as_numpy_iterator():
                y_true.append(y)
            y_true = np.concatenate(y_true, axis=0)
            return y_true
        else:
            return np.array([])

    def evaluate_any(self, idx, batch_size: int = 128, max_steps: int = np.inf, weighted: bool = False):
        """
        Evaluate the custom model on any local data.

        Defaults to run on full data if idx is None.

        :param idx: Indices of observations to evaluate on. Evaluates on all observations if None.
        :param batch_size: Batch size for evaluation.
        :param max_steps: Maximum steps before evaluation round is considered complete.
        :param weighted: Whether to use class weights in evaluation.
        :return: Dictionary of metric names and values.
        """
        if idx is None or idx.any():   # true if the array is not empty or if the passed value is None
            idx = np.arange(0, self.data.n_obs) if idx is None else idx
            dataset = self._get_dataset(
                idx=idx,
                batch_size=batch_size,
                mode='eval',
                weighted=weighted,
                retrieval_batch_size=128,
                shuffle_buffer_size=0,
            )
            results = self.model.training_model.evaluate(x=dataset)
            return dict(zip(self.model.training_model.metrics_names, results))
        else:
            return {}

    def evaluate(self, batch_size: int = 128, max_steps: int = np.inf, weighted: bool = False):
        """
        Evaluate the custom model on local data.

        Defaults to run on full data if idx_test was not set before, ie. train() has not been called before.

        :param batch_size: Batch size for evaluation.
        :param max_steps: Maximum steps before evaluation round is considered complete.
        :param weighted: Whether to use class weights in evaluation.
        :return: Dictionary of metric names and values.
        """
        return self.evaluate_any(idx=self.idx_test, batch_size=batch_size, max_steps=max_steps, weighted=weighted)

    def compute_gradients_input(
            self,
            test_data: bool = False,
            abs_gradients: bool = True
    ):

        if test_data:
            idx = self.idx_test
            n_obs = len(self.idx_test)
        else:
            idx = None
            n_obs = self.data.X.shape[0]

        ds = self._get_dataset(
            idx=idx,
            batch_size=64,
            mode='train_val'  # to get a tf.GradientTape compatible data set
        )
        grads_x = 0
        # Loop over sub-selected data set and sum gradients across all selected observations.
        model = tf.keras.Model(
            self.model.training_model.input,
            self.model.training_model.output
        )

        for step, (x_batch, _, _) in enumerate(ds):
            print("compute gradients wrt. input: batch %i / %i." % (step + 1, np.ceil(n_obs / 64)))
            x = x_batch
            with tf.GradientTape(persistent=True) as tape:
                tape.watch(x)
                model_out = model(x)
            if abs_gradients:
                def f(x):
                    return abs(x)
            else:
                def f(x):
                    return x
            # marginalize on batch level and then accumulate batches
            # batch_jacobian gives output of size: (batch_size, latent_dim, input_dim)
            batch_gradients = f(tape.batch_jacobian(model_out, x).numpy())
            grads_x += np.sum(batch_gradients, axis=0)
        return grads_x
