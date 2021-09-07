import abc
import numpy as np
try:
    import tensorflow as tf
except ImportError:
    tf = None
from sfaira.models.base import BasicModelKeras


class BasicModelEmbedding:

    @abc.abstractmethod
    def predict_reconstructed(self, x, **kwargs):
        pass

    @abc.abstractmethod
    def predict_embedding(self, x, **kwargs):
        pass


class BasicModelKerasEmbedding(BasicModelKeras, BasicModelEmbedding):
    """
    This base class defines model attributes shared across all tf.keras embedding models.
    """

    encoder_model: tf.keras.Model

    def predict_reconstructed(self, x, **kwargs):
        return np.split(self.training_model.predict(x), indices_or_sections=2, axis=1)[0]

    def predict_embedding(self, x, **kwargs):
        return self.encoder_model.predict(x)
