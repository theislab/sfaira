import abc
import anndata
import hashlib
import numpy as np
import pandas
import scipy.sparse
try:
    import tensorflow as tf
except ImportError:
    tf = None
from typing import Union
import os
import warnings
from tqdm import tqdm

from sfaira.consts import AdataIdsSfaira
from sfaira.models import BasicModel
from sfaira.versions.metadata import CelltypeUniverse
from sfaira.versions.topologies import Topologies
from .losses import LossLoglikelihoodNb, LossLoglikelihoodGaussian, LossCrossentropyAgg, KLLoss
from .metrics import custom_mse, custom_negll_nb, custom_negll_gaussian, custom_kl, \
    CustomAccAgg, CustomF1Classwise, CustomFprClasswise, CustomTprClasswise, custom_cce_agg


class EstimatorKeras:
    """
    Estimator base class for keras models.
    """
    data: Union[anndata.AnnData]
    obs_train: Union[pandas.DataFrame, None]
    obs_eval: Union[pandas.DataFrame, None]
    obs_test: Union[pandas.DataFrame, None]
    model: Union[BasicModel, None]
    model_topology: Union[str, None]
    model_id: Union[str, None]
    weights: Union[np.ndarray, None]
    model_dir: Union[str, None]
    history: Union[dict, None]
    train_hyperparam: Union[dict, None]
    idx_train: Union[np.ndarray, None]
    idx_eval: Union[np.ndarray, None]
    idx_test: Union[np.ndarray, None]

    def __init__(
            self,
            data: Union[anndata.AnnData, np.ndarray],
            model_dir: Union[str, None],
            model_id: Union[str, None],
            model_class: Union[str, None],
            organism: Union[str, None],
            organ: Union[str, None],
            model_type: Union[str, None],
            model_topology: Union[str, None],
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', '')
    ):
        self.data = data
        self.obs_train = None
        self.obs_eval = None
        self.obs_test = None
        self.model = None
        self.model_dir = model_dir
        self.model_id = model_id
        self.model_class = model_class.lower()
        self.organism = organism.lower()
        self.organ = organ.lower()
        self.model_type = model_type.lower()
        self.model_topology = model_topology
        self.topology_container = Topologies(
            organism=organism,
            model_class=model_class,
            model_type=model_type,
            topology_id=model_topology
        )

        self.history = None
        self.train_hyperparam = None
        self.idx_train = None
        self.idx_eval = None
        self.idx_test = None
        self.md5 = weights_md5
        self.cache_path = cache_path
        self._adata_ids = AdataIdsSfaira()

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
                        urllib.request.urlretrieve(urljoin(self.model_dir, f'{self.model_id}_weights.data-00000-of-00001'),
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

    def _assert_md5_sum(
            self,
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
            prefetch: int,
            weighted: bool
    ):
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

    def _prepare_data_matrix(self, idx: Union[np.ndarray, None]):
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
            if 'mapped_features' in self.data.uns_keys():
                if self.data.uns[self._adata_ids.mapped_features] == self.topology_container.genome_container.assembly:
                    print(f"found {x.shape[0]} observations")
                    return x

            # Compute indices of genes to keep
            data_ids = self.data.var[self._adata_ids.gene_id_ensembl].values
            idx_feature_kept = np.where([x in self.topology_container.genome_container.ensembl for x in data_ids])[0]
            idx_feature_map = np.array([self.topology_container.genome_container.ensembl.index(x)
                                        for x in data_ids[idx_feature_kept]])

            # Convert to csc and remove unmapped genes
            x = x.tocsc()
            x = x[:, idx_feature_kept]

            # Create reordered feature matrix based on reference and convert to csr
            x_new = scipy.sparse.csc_matrix((x.shape[0], self.topology_container.ngenes), dtype=x.dtype)
            # copying this over to the new matrix in chunks of size `steps` prevents a strange scipy error:
            # ... scipy/sparse/compressed.py", line 922, in _zero_many i, j, offsets)
            # ValueError: could not convert integer scalar
            step = 2000
            if step < len(idx_feature_map):
                for i in range(0, len(idx_feature_map), step):
                    x_new[:, idx_feature_map[i:i + step]] = x[:, i:i + step]
                x_new[:, idx_feature_map[i + step:]] = x[:, i + step:]
            else:
                x_new[:, idx_feature_map] = x

            x_new = x_new.tocsr()

            print(f"found {len(idx_feature_kept)} intersecting features between {x.shape[1]} "
                  f"features in input data set and {self.topology_container.ngenes} features in reference genome")
            print(f"found {x_new.shape[0]} observations")

            return x_new

    def _prepare_sf(self, x):
        if len(x.shape) == 2:
            sf = np.asarray(x.sum(axis=1)).flatten()
        elif len(x.shape) == 1:
            sf = np.asarray(x.sum()).flatten()
        else:
            raise ValueError("x.shape > 2")
        sf = np.log(sf / 1e4 + 1e-10)
        return sf

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
            shuffle_buffer_size: int = int(1e4),
            log_dir: Union[str, None] = None,
            callbacks: Union[list, None] = None,
            weighted: bool = True,
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
        cbs = []
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

        # Split training and evaluation data.
        np.random.seed(1)
        all_idx = np.arange(0, self.data.shape[0])
        if isinstance(test_split, float) or isinstance(test_split, int):
            self.idx_test = np.random.choice(
                a=all_idx,
                size=round(self.data.shape[0] * test_split),
                replace=False,
            )
        elif isinstance(test_split, dict):
            in_test = np.ones((self.data.obs.shape[0],), dtype=int) == 1
            for k, v in test_split.items():
                if isinstance(v, list):
                    in_test = np.logical_and(in_test, np.array([x in v for x in self.data.obs[k].values]))
                else:
                    in_test = np.logical_and(in_test, self.data.obs[k].values == v)
            self.idx_test = np.where(in_test)[0]
            print(f"Found {len(self.idx_test)} out of {self.data.n_obs} cells that correspond to held out data set")
            print(self.idx_test)
        else:
            raise ValueError("type of test_split %s not recognized" % type(test_split))
        idx_train_eval = np.array([x for x in all_idx if x not in self.idx_test])
        np.random.seed(1)
        self.idx_eval = np.random.choice(
            a=idx_train_eval,
            size=round(len(idx_train_eval) * validation_split),
            replace=False
        )
        self.idx_train = np.array([x for x in idx_train_eval if x not in self.idx_eval])

        # Check that none of the train, test, eval partitions are empty
        if not len(self.idx_test):
            warnings.warn("Test partition is empty!")
        if not len(self.idx_eval):
            raise ValueError("The evaluation dataset is empty.")
        if not len(self.idx_train):
            raise ValueError("The train dataset is empty.")

        self.obs_train = self.data.obs.iloc[self.idx_train, :].copy()
        self.obs_eval = self.data.obs.iloc[self.idx_eval, :].copy()
        self.obs_test = self.data.obs.iloc[self.idx_test, :].copy()

        self._compile_models(optimizer=optim)
        train_dataset = self._get_dataset(
            idx=self.idx_train,
            batch_size=batch_size,
            mode='train',
            shuffle_buffer_size=min(shuffle_buffer_size, len(self.idx_train)),
            weighted=weighted
        )
        eval_dataset = self._get_dataset(
            idx=self.idx_eval,
            batch_size=validation_batch_size,
            mode='train_val',
            shuffle_buffer_size=min(shuffle_buffer_size, len(self.idx_eval)),
            weighted=weighted
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

    def get_citations(self):
        """
        Return papers to cite when using this model.

        :return:
        """
        raise NotImplementedError()


class EstimatorKerasEmbedding(EstimatorKeras):
    """
    Estimator class for the embedding model.
    """

    def __init__(
            self,
            data: Union[anndata.AnnData, np.ndarray],
            model_dir: Union[str, None],
            model_id: Union[str, None],
            organism: Union[str, None],
            organ: Union[str, None],
            model_type: Union[str, None],
            model_topology: Union[str, None],
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', '')
    ):
        super(EstimatorKerasEmbedding, self).__init__(
            data=data,
            model_dir=model_dir,
            model_id=model_id,
            model_class="embedding",
            organism=organism,
            organ=organ,
            model_type=model_type,
            model_topology=model_topology,
            weights_md5=weights_md5,
            cache_path=cache_path
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
            raise ValueError('unknown model type %s for EstimatorKerasEmbedding' % self.model_type)
        self.model = Model(
            topology_container=self.topology_container,
            override_hyperpar=override_hyperpar
        )

    @staticmethod
    def _get_output_dim(n_features, model_type, mode='train'):
        if mode == 'predict':  # Output shape is same for predict mode regardless of model type
            output_types = (tf.float32, tf.float32)
            output_shapes = (n_features, ())
        elif model_type == "vae":
            output_types = ((tf.float32, tf.float32), (tf.float32, tf.float32))
            output_shapes = ((n_features, ()), (n_features, ()))
        else:
            output_types = ((tf.float32, tf.float32), tf.float32)
            output_shapes = ((n_features, ()), n_features)

        return output_types, output_shapes

    def _get_dataset(
            self,
            idx: Union[np.ndarray, None],
            batch_size: Union[int, None],
            mode: str,
            shuffle_buffer_size: int = int(1e7),
            prefetch: int = 10,
            weighted: bool = False,
    ):
        """

        :param idx:
        :param batch_size:
        :param mode:
        :param shuffle_buffer_size:
        :param weighted: Whether to use weights. Not implemented for embedding models yet.
        :return:
        """
        # Determine model type [ae, vae(iaf, vamp)]
        model_type = "vae" if self.model_type[:3] == "vae" else "ae"

        if idx is None:
            idx = np.arange(0, self.data.n_obs)

        if mode in ['train', 'train_val', 'eval', 'predict']:
            # Prepare data reading according to whether anndata is backed or not:
            x = self.data.X if self.data.isbacked else self._prepare_data_matrix(idx=idx)

            def generator():
                is_sparse = isinstance(x[0, :], scipy.sparse.spmatrix)
                indices = idx if self.data.isbacked else range(x.shape[0])
                for i in indices:
                    x_sample = x[i, :].toarray().flatten() if is_sparse else x[i, :].flatten()
                    sf = self._prepare_sf(x=x_sample)[0]
                    if mode == 'predict':  # If predicting, only return X regardless of model type
                        yield x_sample, sf
                    elif model_type == "vae":
                        yield (x_sample, sf), (x_sample, sf)
                    else:
                        yield (x_sample, sf), x_sample

            n_features = x.shape[1]
            n_samples = x.shape[0]
            output_types, output_shapes = self._get_output_dim(n_features, model_type, mode=mode)

            dataset = tf.data.Dataset.from_generator(
                generator=generator,
                output_types=output_types,
                output_shapes=output_shapes
            )
            # Only shuffle in train modes
            if mode in ['train', 'train_val']:
                dataset = dataset.repeat()
                dataset = dataset.shuffle(
                    buffer_size=min(n_samples, shuffle_buffer_size),
                    seed=None,
                    reshuffle_each_iteration=True)
            dataset = dataset.batch(batch_size).prefetch(prefetch)

            return dataset

        elif mode == 'gradient_method':
            # Prepare data reading according to whether anndata is backed or not:
            if self.data.isbacked:
                n_features = self.data.X.shape[1]
                cell_to_class = self._get_class_dict(obs_key=self._adata_ids.cell_ontology_class)
                output_types, output_shapes = self._get_output_dim(n_features, 'vae')

                def generator():
                    sparse = isinstance(self.data.X[0, :], scipy.sparse.spmatrix)
                    for i in idx:
                        x = self.data.X[i, :].toarray().flatten() if sparse else self.data.X[i, :].flatten()
                        sf = self._prepare_sf(x=x)[0]
                        y = self.data.obs[self._adata_ids.cell_ontology_class][i]
                        yield (x, sf), (x, cell_to_class[y])
            else:
                x = self._prepare_data_matrix(idx=idx)
                sf = self._prepare_sf(x=x)
                cell_to_class = self._get_class_dict(obs_key=self._adata_ids.cell_ontology_class)
                y = self.data.obs[self._adata_ids.cell_ontology_class][idx]  # for gradients per celltype in compute_gradients_input()
                n_features = x.shape[1]
                output_types, output_shapes = self._get_output_dim(n_features, 'vae')

                def generator():
                    for i in range(x.shape[0]):
                        yield (x[i, :].toarray().flatten(), sf[i]), (x[i, :].toarray().flatten(), cell_to_class[y[i]])

            dataset = tf.data.Dataset.from_generator(
                generator=generator,
                output_types=output_types,
                output_shapes=output_shapes
            )
            dataset = dataset.shuffle(
                buffer_size=shuffle_buffer_size,
                seed=None,
                reshuffle_each_iteration=True
            ).batch(batch_size).prefetch(prefetch)

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

    def evaluate_any(self, idx, batch_size=64, max_steps=20):
        """
        Evaluate the custom model on any local data.

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
                mode='eval'
            )
            steps = min(max(len(idx) // batch_size, 1), max_steps)
            results = self.model.training_model.evaluate(
                x=dataset, steps=steps
            )
            return dict(zip(self.model.training_model.metrics_names, results))
        else:
            return {}

    def evaluate(self, batch_size=64, max_steps=20):
        """
        Evaluate the custom model on local data.

        Defaults to run on full data if idx_test was not set before, ie. train() has not been called before.

        :return: Dictionary of metric names and values.
        """
        if self.idx_test is None or self.idx_test.any():  # true if the array is not empty or if the passed value is None
            idx = np.arange(0, self.data.n_obs) if self.idx_test is None else self.idx_test
            dataset = self._get_dataset(
                idx=idx,
                batch_size=batch_size,
                mode='eval'
            )
            steps = min(max(len(idx) // batch_size, 1), max_steps)
            results = self.model.training_model.evaluate(
                x=dataset, steps=steps
            )
            return dict(zip(self.model.training_model.metrics_names, results))
        else:
            return {}

    def predict(self):
        """
        return the prediction of the model

        :return:
        prediction
        """
        if self.idx_test is None or self.idx_test.any():  # true if the array is not empty or if the passed value is None
            x = self._get_dataset(
                idx=self.idx_test,
                batch_size=64,
                mode='predict'
            )
            return self.model.predict_reconstructed(
                x=x
            )
        else:
            return np.array([])

    def predict_embedding(self):
        """
        return the prediction in the latent space (z_mean for variational models)

        :return:
        latent space
        """
        if self.idx_test is None or self.idx_test.any():  # true if the array is not empty or if the passed value is None
            x = self._get_dataset(
                idx=self.idx_test,
                batch_size=64,
                mode='predict'
            )
            return self.model.predict_embedding(
                x=x,
                variational=False
            )
        else:
            return np.array([])

    def predict_embedding_variational(self):
        """
        return the prediction of z, z_mean, z_log_var in the variational latent space

        :return:
        sample of latent space, mean of latent space, variance of latent space
        """
        if self.idx_test is None or self.idx_test:  # true if the array is not empty or if the passed value is None
            x = self._get_dataset(
                idx=self.idx_test,
                batch_size=64,
                mode='predict'
            )
            return self.model.predict_embedding(
                x=x,
                variational=True
            )
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
            cell_to_id = self._get_class_dict(obs_key=self._adata_ids.cell_ontology_class)
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
                def f(x):
                    return abs(x)
            else:
                def f(x):
                    return x
            # marginalize on batch level and then accumulate batches
            # batch_jacobian gives output of size: (batch_size, latent_dim, input_dim)
            batch_gradients = f(tape.batch_jacobian(model_out, x))
            return batch_gradients

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

    celltypes_version: CelltypeUniverse

    def __init__(
            self,
            data: Union[anndata.AnnData, np.ndarray],
            model_dir: Union[str, None],
            model_id: Union[str, None],
            organism: Union[str, None],
            organ: Union[str, None],
            model_type: Union[str, None],
            model_topology: Union[str, None],
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', ''),
            max_class_weight: float = 1e3
    ):
        super(EstimatorKerasCelltype, self).__init__(
            data=data,
            model_dir=model_dir,
            model_id=model_id,
            model_class="celltype",
            organism=organism,
            organ=organ,
            model_type=model_type,
            model_topology=model_topology,
            weights_md5=weights_md5,
            cache_path=cache_path
        )
        self.max_class_weight = max_class_weight
        self.celltypes_version = CelltypeUniverse(organism=organism)

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
            celltypes_version=self.celltypes_version,
            topology_container=self.topology_container,
            override_hyperpar=override_hyperpar
        )

    @property
    def ids(self):
        return self.celltypes_version.target_universe

    @property
    def ntypes(self):
        return self.celltypes_version.ntypes

    @property
    def ontology_ids(self):
        return self.celltypes_version.target_universe

    def _get_celltype_out(
            self,
            idx: Union[np.ndarray, None],
            lookup_ontology="names"
    ):
        """
        Build one hot encoded cell type output tensor and observation-wise weight matrix.

        :param lookup_ontology: list of ontology names to consider.
        :return:
        """
        if idx is None:
            idx = np.arange(0, self.data.n_obs)
        # One whether "unknown" is already included, otherwise add one extra column.
        if np.any([x.lower() == "unknown" for x in self.ids]):
            type_classes = self.ntypes
        else:
            type_classes = self.ntypes + 1
        y = np.zeros((len(idx), type_classes), dtype="float32")
        celltype_idx = self.model.celltypes_version.map_to_target_leaves(
            nodes=self.data.obs[self._adata_ids.cell_ontology_class].values[idx].tolist(),
            ontology="custom",
            ontology_id=lookup_ontology,
            return_type="idx"
        )
        for i, x in enumerate(celltype_idx):
            # Distribute probability mass uniformly across classes if multiple classes match:
            y[i, x] = 1. / len(x)
        # Distribute aggregated class weight for computation of weights:
        freq = np.mean(y / np.sum(y, axis=1, keepdims=True), axis=0, keepdims=True)
        weights = 1. / np.matmul(y, freq.T)  # observation wise weight matrix
        # Threshold weights:
        weights = np.asarray(
            np.minimum(weights, np.zeros_like(weights) + self.max_class_weight),
            dtype="float32"
        ).flatten()
        return weights, y

    def _get_dataset(
            self,
            idx: Union[np.ndarray, None],
            batch_size: Union[int, None],
            mode: str,
            shuffle_buffer_size: int = int(1e7),
            prefetch: int = 10,
            weighted: bool = True,
    ):
        """

        :param idx:
        :param batch_size:
        :param mode:
        :param shuffle_buffer_size:
        :param weighted: Whether to use weights.
        :return:
        """
        if mode == 'train' or mode == 'train_val':
            weights, y = self._get_celltype_out(idx=idx)
            if not weighted:
                weights = np.ones_like(weights)

            if self.data.isbacked:
                n_features = self.data.X.shape[1]

                def generator():
                    sparse = isinstance(self.data.X[0, :], scipy.sparse.spmatrix)
                    for i, ii in enumerate(idx):
                        x = self.data.X[ii, :].toarray().flatten() if sparse else self.data.X[ii, :].flatten()
                        yield x, y[i, :], weights[i]
            else:
                x = self._prepare_data_matrix(idx=idx)
                n_features = x.shape[1]

                def generator():
                    for i, ii in enumerate(idx):
                        yield x[i, :].toarray().flatten(), y[i, :], weights[i]

            dataset = tf.data.Dataset.from_generator(
                generator=generator,
                output_types=(tf.float32, tf.float32, tf.float32),
                output_shapes=(
                    (tf.TensorShape([n_features])),
                    tf.TensorShape([y.shape[1]]),
                    tf.TensorShape([])
                )
            )
            if mode == 'train':
                dataset = dataset.repeat()
            dataset = dataset.shuffle(
                buffer_size=min(x.shape[0], shuffle_buffer_size),
                seed=None,
                reshuffle_each_iteration=True
            ).batch(batch_size).prefetch(prefetch)

            return dataset

        elif mode == 'eval':
            weights, y = self._get_celltype_out(idx=idx)
            if not weighted:
                weights = np.ones_like(weights)

            # Prepare data reading according to whether anndata is backed or not:
            if self.data.isbacked:
                # Need to supply sorted indices to backed anndata:
                x = self.data.X[np.sort(idx), :]
                # Sort back in original order of indices.
                x = x[[np.where(np.sort(idx) == i)[0][0] for i in idx], :]
            else:
                x = self._prepare_data_matrix(idx=idx)
                x = x.toarray()

            return x, y, weights

        elif mode == 'predict':
            # Prepare data reading according to whether anndata is backed or not:
            if self.data.isbacked:
                # Need to supply sorted indices to backed anndata:
                x = self.data.X[np.sort(idx), :]
                # Sort back in original order of indices.
                x = x[[np.where(np.sort(idx) == i)[0][0] for i in idx], :]
            else:
                x = self._prepare_data_matrix(idx=idx)
                x = x.toarray()

            return x, None, None

        else:
            raise ValueError(f'Mode {mode} not recognised. Should be "train", "eval" or" predict"')

    def _get_loss(self):
        return LossCrossentropyAgg()

    def _metrics(self):
        if np.any([x.lower() == "unknown" for x in self.ids]):
            ntypes = self.ntypes
        else:
            ntypes = self.ntypes + 1
        return [
            "accuracy",
            custom_cce_agg,
            CustomAccAgg(),
            CustomF1Classwise(k=ntypes),
            CustomFprClasswise(k=ntypes),
            CustomTprClasswise(k=ntypes)
        ]

    def predict(self):
        """
        Return the prediction of the model

        :return:
        prediction
        """
        if self.idx_test is None or self.idx_test.any():   # true if the array is not empty or if the passed value is None
            x, _, _ = self._get_dataset(
                idx=self.idx_test,
                batch_size=None,
                mode='predict'
            )
            return self.model.training_model.predict(
                x=x
            )
        else:
            return np.array([])

    def ytrue(self):
        """
        Return the true labels of the test set.

        :return: true labels
        """
        if self.idx_test is None or self.idx_test.any():   # true if the array is not empty or if the passed value is None
            x, y, w = self._get_dataset(
                idx=self.idx_test,
                batch_size=None,
                mode='eval'
            )
            return y
        else:
            return np.array([])

    def evaluate_any(
            self,
            idx,
            weighted: bool = True
    ):
        """
        Evaluate the custom model on any local data.

        :param idx: Indices of observations to evaluate on. Evaluates on all observations if None.
        :param weighted: Whether to use class weights in evaluation.
        :return: Dictionary of metric names and values.
        """
        if idx is None or idx.any():   # true if the array is not empty or if the passed value is None
            x, y, w = self._get_dataset(
                idx=idx,
                batch_size=None,
                mode='eval',
                weighted=weighted
            )
            results = self.model.training_model.evaluate(
                x=x, y=y, sample_weight=w
            )
            return dict(zip(self.model.training_model.metrics_names, results))
        else:
            return {}

    def evaluate(self, weighted: bool = True):
        """
        Evaluate the custom model on local data.

        Defaults to run on full data if idx_test was not set before, ie. train() has not been called before.

        :param weighted: Whether to use class weights in evaluation.
        :return: model.evaluate
        """
        if self.idx_test is None or self.idx_test.any():   # true if the array is not empty or if the passed value is None
            x, y, w = self._get_dataset(
                idx=self.idx_test,
                batch_size=None,
                mode='eval',
                weighted=weighted
            )
            return self.model.training_model.evaluate(
                x=x, y=y, sample_weight=w
            )
        else:
            return np.array([])

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
