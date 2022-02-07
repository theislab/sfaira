import abc
import hashlib
import os
from abc import ABC
from typing import List, Union

import anndata
import numpy as np
import tensorflow as tf
from tqdm import tqdm

from sfaira.consts import AdataIdsSfaira, AdataIds
from sfaira.data.store.carts.single import CartSingle
from sfaira.data.store.stores.single import StoreSingleFeatureSpace
from sfaira.estimators.base import EstimatorBaseCelltype, EstimatorBaseEmbedding
from sfaira.models import BasicModelKeras
from sfaira.models.celltype import BasicModelKerasCelltype
from sfaira.models.embedding import BasicModelKerasEmbedding
from sfaira.versions.metadata import CelltypeUniverse, OntologyObo
from sfaira.versions.topologies import TopologyContainer
from .losses import LossLoglikelihoodNb, LossLoglikelihoodGaussian, LossCrossentropyAgg, KLLoss
from .metrics import custom_mse, custom_negll_nb, custom_negll_gaussian, \
    CustomAccAgg, CustomF1Classwise, CustomFprClasswise, CustomTprClasswise, custom_cce_agg


def prepare_sf(x):
    """
    Uses a minimal size factor of 1e-3 for total counts / 1e4
    """
    sf = np.asarray(x.sum(axis=1, keepdims=True))
    sf = np.log(np.maximum(sf / 1e4, 1e-3))
    return sf


def process_tf_dataset(dataset, mode: str, batch_size: int, cache: bool, shuffle_buffer_size: int, prefetch)\
        -> tf.data.Dataset:
    if cache:
        dataset = dataset.cache()
    if mode in ['train', 'train_val']:
        dataset = dataset.repeat()
        if shuffle_buffer_size > 0:
            # Only shuffle in train modes
            dataset = dataset.shuffle(buffer_size=shuffle_buffer_size, seed=None, reshuffle_each_iteration=True)
    dataset = dataset.batch(batch_size).prefetch(prefetch)
    return dataset


def get_optimizer(optimizer: str, lr: float):
    if optimizer.lower() == "adam":
        return tf.keras.optimizers.Adam(learning_rate=lr)
    elif optimizer.lower() == "sgd":
        return tf.keras.optimizers.SGD(learning_rate=lr)
    elif optimizer.lower() == "rmsprop":
        return tf.keras.optimizers.RMSprop(learning_rate=lr)
    elif optimizer.lower() == "adagrad":
        return tf.keras.optimizers.Adagrad(learning_rate=lr)
    else:
        assert False


def assemble_cbs(patience, lr_schedule_factor, lr_schedule_patience, lr_schedule_min_lr, lr_verbose, log_dir,
                 callbacks) -> List[tf.keras.callbacks.Callback]:
    cbs = [tf.keras.callbacks.TerminateOnNaN()]
    if patience is not None and patience > 0:
        cbs.append(tf.keras.callbacks.EarlyStopping(
            monitor='val_loss',
            patience=patience,
            restore_best_weights=True,
            verbose=1
        ))
    if lr_schedule_factor is not None and lr_schedule_factor < 1.:
        cbs.append(tf.keras.callbacks.ReduceLROnPlateau(
            monitor='val_loss',
            factor=lr_schedule_factor,
            patience=lr_schedule_patience,
            min_lr=lr_schedule_min_lr,
            verbose=lr_verbose
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
    return cbs


class EstimatorKeras(ABC):
    """
    Estimator base class for keras models.

    Important: Subclass implementing abstract classes also has to inherit from EstimatorBase class.
    """
    data: StoreSingleFeatureSpace
    model: Union[BasicModelKeras, None]
    weights: Union[np.ndarray, None]
    model_dir: Union[str, None]
    history: Union[dict, None]
    train_hyperparam: Union[dict, None]
    idx_train: Union[np.ndarray, None]
    idx_eval: Union[np.ndarray, None]
    idx_test: Union[np.ndarray, None]
    cache_path: str
    model_dir: Union[str, None]
    model_id: Union[str, None]
    md5: Union[str, None]

    @abc.abstractmethod
    def _get_cart(self, **kwargs) -> CartSingle:
        """
        Yield a generator based on which a tf dataset can be built.
        """
        pass

    @abc.abstractmethod
    def _tf_dataset_kwargs(self, mode: str):
        pass

    @abc.abstractmethod
    def split_train_val_test(self, val_split: float, test_split: Union[float, dict]):
        """
        Split indices in store into train, valiation and test split.
        """
        pass

    def get_one_time_tf_dataset(self, idx, mode, batch_size=None, prefetch=None):
        batch_size = 128 if batch_size is None else batch_size
        prefetch = 10 if prefetch is None else prefetch
        tf_kwargs = {
            "batch_size": batch_size,
            "cache": False,
            "prefetch": prefetch,
            "shuffle_buffer_size": 0,
        }
        train_gen = self._get_cart(idx=idx, mode=mode, retrieval_batch_size=128,
                                   randomized_batch_access=False)
        train_tf_dataset_kwargs = self._tf_dataset_kwargs(mode=mode)
        train_dataset = train_gen.adaptor(generator_type="tensorflow", **train_tf_dataset_kwargs)
        train_dataset = process_tf_dataset(dataset=train_dataset, mode=mode, **tf_kwargs)
        return train_dataset

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
        # Save training settings to allow model restoring.
        self.train_hyperparam = {"epochs": epochs, "max_steps_per_epoch": max_steps_per_epoch, "optimizer": optimizer,
                                 "lr": lr, "batch_size": batch_size, "validation_split": validation_split,
                                 "validation_batch_size": validation_batch_size,
                                 "max_validation_steps": max_validation_steps, "patience": patience,
                                 "lr_schedule_min_lr": lr_schedule_min_lr, "lr_schedule_factor": lr_schedule_factor,
                                 "lr_schedule_patience": lr_schedule_patience, "log_dir": log_dir, "weighted": weighted}

        # Set callbacks.
        cbs = assemble_cbs(patience=patience, lr_schedule_factor=lr_schedule_factor,
                           lr_schedule_patience=lr_schedule_patience, lr_schedule_min_lr=lr_schedule_min_lr,
                           lr_verbose=verbose, log_dir=log_dir, callbacks=callbacks)

        # Check randomisation settings:
        if shuffle_buffer_size is not None and shuffle_buffer_size > 0 and randomized_batch_access:
            raise ValueError("You are using shuffle_buffer_size and randomized_batch_access, this is likely not "
                             "intended.")
        shuffle_buffer_size = shuffle_buffer_size if shuffle_buffer_size is not None else 0
        if cache_full and randomized_batch_access:
            raise ValueError("You are using cache_full and randomized_batch_access, this is likely not intended.")
        self.split_train_val_test(val_split=validation_split, test_split=test_split)
        self._compile_models(optimizer=get_optimizer(optimizer=optimizer, lr=lr))

        tf_kwargs = {"batch_size": batch_size, "cache": cache_full, "prefetch": prefetch,
                     "shuffle_buffer_size": min(shuffle_buffer_size, len(self.idx_train))}
        train_gen = self._get_cart(idx=self.idx_train, mode='train', retrieval_batch_size=retrieval_batch_size,
                                   randomized_batch_access=randomized_batch_access, weighted=weighted)
        train_tf_dataset_kwargs = self._tf_dataset_kwargs(mode="train")
        train_dataset = train_gen.adaptor(generator_type="tensorflow", **train_tf_dataset_kwargs)
        train_dataset = process_tf_dataset(dataset=train_dataset, mode="train", **tf_kwargs)
        val_gen = self._get_cart(idx=self.idx_train, mode='train', retrieval_batch_size=retrieval_batch_size,
                                 randomized_batch_access=randomized_batch_access, weighted=weighted)
        val_tf_dataset_kwargs = self._tf_dataset_kwargs(mode="train_val")
        val_dataset = val_gen.adaptor(generator_type="tensorflow", **val_tf_dataset_kwargs)
        val_dataset = process_tf_dataset(dataset=val_dataset, mode="train", **tf_kwargs)

        steps_per_epoch = min(max(len(self.idx_train) // batch_size, 1), max_steps_per_epoch)
        validation_steps = min(max(len(self.idx_eval) // validation_batch_size, 1), max_validation_steps)

        self.history = self.model.training_model.fit(
            x=train_dataset,
            epochs=epochs,
            steps_per_epoch=steps_per_epoch,
            callbacks=cbs,
            validation_data=val_dataset,
            validation_steps=validation_steps,
            verbose=verbose
        ).history

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
                                           os.path.join(self.cache_path, os.path.basename(self.model_dir)))
                fn = os.path.join(self.cache_path, os.path.basename(self.model_dir))
            except HTTPError:
                try:
                    urllib.request.urlretrieve(urljoin(self.model_dir, f'{self.model_id}_weights.h5'),
                                               os.path.join(self.cache_path, f'{self.model_id}_weights.h5'))
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

        if self.md5 is not None:
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

    @staticmethod
    def _assert_md5_sum(
            fn,
            target_md5
    ):
        with open(fn, 'rb') as f:
            hsh = hashlib.md5(f.read()).hexdigest()
        if not hsh == target_md5:
            raise ValueError("md5 of %s did not match expectation" % fn)

    def init_model(self, clear_weight_cache=True, override_hyperpar=None):
        """
        Instantiate the model.
        :return:
        """
        if clear_weight_cache:
            if os.path.exists(os.path.join(self.cache_path, 'weights')):
                for file in os.listdir(os.path.join(self.cache_path, 'weights')):
                    file_path = os.path.join(os.path.join(self.cache_path, 'weights'), file)
                    os.remove(file_path)


class EstimatorKerasEmbedding(EstimatorBaseEmbedding, EstimatorKeras):
    """
    Estimator class for the embedding model.
    """

    model: Union[BasicModelKerasEmbedding, None]

    def __init__(
            self,
            data: Union[anndata.AnnData, np.ndarray, StoreSingleFeatureSpace],
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

    def _tf_dataset_kwargs(self, mode: str):
        # Determine model type [ae, vae(iaf, vamp)]
        model_type = "vae" if self.model_type[:3] == "vae" else "ae"
        output_types_x = (tf.float32, tf.float32)
        output_shapes_x = (self.data.n_vars, 1)
        if mode == 'predict':  # Output shape is same for predict mode regardless of model type
            output_types = output_types_x,
            output_shapes = output_shapes_x,
        elif model_type == "vae":
            output_types = (output_types_x, (tf.float32, tf.float32))
            output_shapes = (output_shapes_x, (self.data.n_vars, 1))
        else:
            output_types = (output_types_x, (tf.float32, ))
            output_shapes = (output_shapes_x, (self.data.n_vars, ))
        return {"output_types": output_types, "output_shapes": output_shapes}

    def _get_cart(
            self,
            idx: Union[np.ndarray, None],
            mode: str,
            retrieval_batch_size: int,
            randomized_batch_access: bool,
            **kwargs
    ):
        # Define constants used by map_fn in outer name space so that they are not created for each sample.
        model_type = "vae" if self.model_type[:3] == "vae" else "ae"

        def map_fn(x_sample, obs_sample):
            x_sample = np.asarray(x_sample)
            sf_sample = prepare_sf(x=x_sample)
            output_x = (x_sample, sf_sample)
            if mode == 'predict':
                output = output_x,
            elif model_type == "vae":
                output = output_x, (x_sample, sf_sample)
            else:
                output = output_x, (x_sample, )
            return output

        g = self.data.checkout(idx=idx, retrieval_batch_size=retrieval_batch_size, obs_keys=[], map_fn=map_fn,
                               return_dense=True, randomized_batch_access=randomized_batch_access,
                               random_access=False)
        return g

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
        idx = self._process_idx_for_eval(idx=idx)
        if idx is not None:
            dataset = self.get_one_time_tf_dataset(idx=idx, batch_size=batch_size, mode='eval')
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
        idx = self._process_idx_for_eval(idx=self.idx_test)
        if idx is not None:
            return self.evaluate_any(idx=self.idx_test, batch_size=batch_size, max_steps=max_steps)
        else:
            return {}

    def predict(self, batch_size: int = 128):
        """
        return the prediction of the model

        :return:
        prediction
        """
        idx = self._process_idx_for_eval(idx=self.idx_test)
        if idx is not None:
            dataset = self.get_one_time_tf_dataset(idx=idx, batch_size=batch_size, mode='predict')
            return self.model.predict_reconstructed(x=dataset)
        else:
            return np.array([])

    def predict_embedding(self, batch_size: int = 128, variational: bool = False):
        """
        return the prediction in the latent space (z_mean for variational models)

        :params variational: Whether toreturn the prediction of z, z_mean, z_log_var in the variational latent space.
        :return:
        latent space
        """
        idx = self._process_idx_for_eval(idx=self.idx_test)
        if len(idx) > 0:
            dataset = self.get_one_time_tf_dataset(idx=idx, batch_size=batch_size, mode='predict')
            return self.model.predict_embedding(x=dataset, variational=variational)
        else:
            return np.array([])

    def compute_gradients_input(
            self,
            batch_size: int = 128,
            test_data: bool = False,
            abs_gradients: bool = True,
            per_celltype: bool = False
    ):
        # TODO may need to be adapted to new dataset / generator format
        if test_data:
            idx = self.idx_test
            if self.idx_test is None:
                num_samples = 10000
                idx = np.random.randint(0, self.data.n_obs, num_samples)
            n_obs = len(idx)
        else:
            idx = None
            n_obs = self.data.n_obs

        ds = self.get_one_time_tf_dataset(idx=idx, batch_size=batch_size, mode='gradient_method')

        if per_celltype:
            cell_to_id = self._get_class_dict(obs_key=self._adata_ids.cell_type)
            cell_names = cell_to_id.keys()
            cell_id = cell_to_id.values()
            id_to_cell = dict([(key, value) for (key, value) in zip(cell_id, cell_names)])
            grads_x = dict([(key, 0) for key in cell_names])
            counts = dict([(key, 0) for key in cell_names])
        else:
            grads_x = 0
        # Loop over sub-selected data set and sum gradients across all selected observations.
        # TODO: fix bug for vaeiaf model. This function can not be called for vaeiaf model
        if self.model_type[:3] == "vae":
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


class EstimatorKerasCelltype(EstimatorBaseCelltype, EstimatorKeras):
    """
    Estimator class for the cell type model.
    """

    celltype_universe: CelltypeUniverse
    model: Union[BasicModelKerasCelltype, None]

    def __init__(
            self,
            data: Union[anndata.AnnData, StoreSingleFeatureSpace],
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
            model_id=model_id,
            model_topology=model_topology,
            weights_md5=weights_md5,
            cache_path=cache_path,
            celltype_ontology=celltype_ontology,
            max_class_weight=max_class_weight,
            remove_unlabeled_cells=remove_unlabeled_cells,
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

    def _tf_dataset_kwargs(self, mode):
        output_types_x = (tf.float32,)
        output_shapes_x = (tf.TensorShape([self.data.n_vars]), )
        if mode == 'predict':
            output_types = (output_types_x, )
            output_shapes = (output_shapes_x, )
        else:
            output_types = (output_types_x, tf.float32)
            output_shapes = (output_shapes_x, tf.TensorShape([self.ntypes]))
        return {"output_types": output_types, "output_shapes": output_shapes}

    def _get_cart(
            self,
            idx: Union[np.ndarray, None],
            mode: str,
            retrieval_batch_size: int,
            randomized_batch_access: bool,
            weighted: bool = False,
            **kwargs
    ) -> CartSingle:
        # Define constants used by map_fn in outer name space so that they are not created for each sample.
        if weighted:
            raise ValueError("using weights with store is not supported yet")
        yield_labels = mode in ["train", "train_val", "eval", "test"]
        if yield_labels:
            onehot_encoder = self._one_hot_encoder()

        def map_fn(x_sample, obs_sample):
            x_sample = np.asarray(x_sample)
            output_x = (x_sample, )
            if yield_labels:
                y_sample = onehot_encoder(obs_sample[self._adata_ids.cell_type + self._adata_ids.onto_id_suffix].values)
                # Only yield observations with valid label:
                idx_keep = y_sample.sum(axis=1) > 0.
                if not np.all(idx_keep):
                    idx_keep = np.where(idx_keep)[0]
                    output_x = tuple([x[idx_keep, :] for x in output_x])
                    y_sample = y_sample[idx_keep, :]
                output_y = y_sample
                output = output_x, output_y
            else:
                output = output_x,
            return output

        g = self.data.checkout(idx=idx, retrieval_batch_size=retrieval_batch_size,
                               obs_keys=[self._adata_ids.cell_type + self._adata_ids.onto_id_suffix], map_fn=map_fn,
                               return_dense=True, randomized_batch_access=randomized_batch_access,
                               random_access=False)
        return g

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
        idx = self._process_idx_for_eval(idx=self.idx_test)
        if len(idx) > 0:
            dataset = self.get_one_time_tf_dataset(idx=idx, batch_size=batch_size, mode='predict')
            return self.model.training_model.predict(x=dataset)
        else:
            return np.array([])

    def ytrue(self, batch_size: int = 128, max_steps: int = np.inf):
        """
        Return the true labels of the test set.

        :return: true labels
        """
        idx = self._process_idx_for_eval(idx=self.idx_test)
        if len(idx) > 0:
            dataset = self.get_one_time_tf_dataset(idx=idx, batch_size=batch_size, mode='eval')
            y_true = []
            for _, y in dataset.as_numpy_iterator():
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
        idx = self._process_idx_for_eval(idx=idx)
        if len(idx) > 0:
            dataset = self.get_one_time_tf_dataset(idx=idx, batch_size=batch_size, mode='eval')
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
        # TODO may need to be adapted to new dataset / generator format
        if test_data:
            idx = self.idx_test
            n_obs = len(self.idx_test)
        else:
            idx = None
            n_obs = self.data.n_obs

        # to get a tf.GradientTape compatible data set
        ds = self.get_one_time_tf_dataset(idx=idx, batch_size=64, mode='train_val')
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
