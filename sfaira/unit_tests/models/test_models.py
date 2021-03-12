import abc
import numpy as np
import tensorflow as tf
import unittest

from sfaira.estimators.losses import LossLoglikelihoodNb
from sfaira.estimators.metrics import custom_mse

import sfaira.models as models
from sfaira.models.base import BasicModel


class _TestModel:
    model: BasicModel
    data: np.ndarray

    @abc.abstractmethod
    def init_model(self):
        """
        Initialise target model as .model attribute.

        :return:
        """
        pass

    def simulate(self):
        """
        Simulate basic data example used for unit test.

        Sets attribute .data with simulated data.

        :return:
        """
        self.data = np.random.uniform(low=0, high=100, size=(1000, 100)).astype('float32')
        self.sf = np.zeros((1000, 1))


class TestModelAe(unittest.TestCase, _TestModel):

    def init_model(self):
        tf.compat.v1.set_random_seed(0)
        self.dataset = tf.data.Dataset.from_tensor_slices(
            ((self.data, self.sf), self.data)
        )
        self.model = models.embedding.ModelAe(in_dim=self.data.shape[1])

    def compile_models(self):
        self.model.training_model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
            loss=LossLoglikelihoodNb(),
            metrics=[custom_mse]
        )

    def train(self):
        self.model.training_model.fit(
            self.dataset.repeat().batch(128),
            epochs=2, steps_per_epoch=100
        )

    def test_for_fatal(self):
        print(tf.__version__)
        np.random.seed(1)
        self.simulate()
        self.init_model()
        self.compile_models()
        self.train()
        _ = self.model.training_model.evaluate(x=(self.data, self.sf), y=self.data)
        embedding = self.model.predict_embedding(x=(self.data, self.sf))
        assert embedding.shape[0] == self.data.shape[0], embedding.shape
        denoised = self.model.predict_reconstructed(x=(self.data, self.sf))
        assert denoised.shape == self.data.shape, (denoised.shape, self.data.shape)
        return True


class TestModelVae(unittest.TestCase, _TestModel):

    def init_model(self):
        # (_,_), (_,sf) is dummy for kl loss
        self.dataset = tf.data.Dataset.from_tensor_slices(
            ((self.data, self.sf), (self.data, self.sf))
        )
        self.model = models.embedding.ModelVae(in_dim=self.data.shape[1])
        tf.compat.v1.set_random_seed(0)

    def compile_models(self):
        self.model.training_model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
            loss=LossLoglikelihoodNb(),
            metrics=[custom_mse]
        )

    def train(self):
        self.model.training_model.fit(
            self.dataset.repeat().batch(128),
            epochs=2, steps_per_epoch=100
        )

    def test_for_fatal(self):
        print(tf.__version__)
        np.random.seed(1)
        self.simulate()
        self.init_model()
        self.compile_models()
        self.train()
        # (_,_), (_,sf) is dummy for kl loss
        _ = self.model.training_model.evaluate(x=(self.data, self.sf), y=(self.data, self.sf))
        embedding = self.model.predict_embedding(x=(self.data, self.sf))
        assert embedding.shape[0] == self.data.shape[0], embedding.shape
        denoised = self.model.predict_reconstructed(x=(self.data, self.sf))
        assert denoised.shape == self.data.shape, (denoised.shape, self.data.shape)
        return True


class TestModelLinear(unittest.TestCase, _TestModel):

    def init_model(self):
        self.dataset = tf.data.Dataset.from_tensor_slices(
            ((self.data, self.sf), self.data)
        )
        self.model = models.embedding.ModelLinear(
            in_dim=self.data.shape[1],
            output_layer="nb_shared_disp"
        )
        tf.compat.v1.set_random_seed(0)

    def compile_models(self):
        self.model.training_model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
            loss=LossLoglikelihoodNb(),
            metrics=[custom_mse]
        )

    def train(self):
        self.model.training_model.fit(
            self.dataset.repeat().batch(128),
            epochs=2, steps_per_epoch=100
        )

    def test_for_fatal(self):
        print(tf.__version__)
        np.random.seed(1)
        self.simulate()
        self.init_model()
        self.compile_models()
        self.train()
        _ = self.model.training_model.evaluate(x=(self.data, self.sf), y=self.data)
        embedding = self.model.predict_embedding(x=(self.data, self.sf))
        assert embedding.shape[0] == self.data.shape[0], (embedding.shape, self.data.shape)
        denoised = self.model.predict_reconstructed(x=(self.data, self.sf))
        assert denoised.shape == self.data.shape, (denoised.shape, self.data.shape)
        return True


class TestCelltypeMlp(unittest.TestCase, _TestModel):

    def init_model(self):
        self.out_dim = 20
        self.dataset = tf.data.Dataset.from_tensor_slices(
            (self.data, np.ones((self.data.shape[0], self.out_dim)))
        )
        self.model = models.celltype.CellTypeMlp(
            in_dim=self.data.shape[1],
            out_dim=self.out_dim,
            units=[30]
        )

    def compile_model(self):
        self.model.training_model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
            loss=tf.keras.losses.CategoricalCrossentropy(),
            metrics=["accuracy"]
        )

    def train(self):
        train_dataset = self.dataset
        train_dataset = train_dataset.shuffle(buffer_size=1024).batch(64)
        self.model.training_model.fit(train_dataset, epochs=5, steps_per_epoch=2)

    def test_for_fatal(self):
        print(tf.__version__)
        np.random.seed(1)
        self.simulate()
        self.init_model()
        self.compile_model()
        self.train()
        self.model.training_model.evaluate(
            x=self.data, y=np.ones((self.data.shape[0], self.out_dim))
        )
        self.model.training_model.predict(
            x=self.data
        )
        return True


class TestCelltypeMarker(unittest.TestCase, _TestModel):

    def init_model(self):
        self.out_dim = 20
        self.dataset = tf.data.Dataset.from_tensor_slices(
            (self.data, np.ones((self.data.shape[0], self.out_dim)))
        )
        self.model = models.celltype.CellTypeMarker(
            in_dim=self.data.shape[1],
            out_dim=self.out_dim,
        )

    def compile_model(self):
        self.model.training_model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
            loss=tf.keras.losses.CategoricalCrossentropy(),
            metrics=["accuracy"]
        )

    def train(self):
        train_dataset = self.dataset
        train_dataset = train_dataset.shuffle(buffer_size=1024).batch(64)
        self.model.training_model.fit(train_dataset, epochs=5, steps_per_epoch=2)

    def test_for_fatal(self):
        print(tf.__version__)
        np.random.seed(1)
        self.simulate()
        self.init_model()
        self.compile_model()
        self.train()
        self.model.training_model.evaluate(
            x=self.data, y=np.ones((self.data.shape[0], self.out_dim))
        )
        self.model.training_model.predict(
            x=self.data
        )
        return True


if __name__ == '__main__':
    unittest.main()
