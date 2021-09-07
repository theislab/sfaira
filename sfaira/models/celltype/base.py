import abc
try:
    import tensorflow as tf
except ImportError:
    tf = None
from sfaira.models.base import BasicModelKeras


class BasicModelEmbedding:

    @abc.abstractmethod
    def predict(self, x, **kwargs):
        pass


class BasicModelKerasCelltype(BasicModelKeras):
    """
    This base class defines model attributes shared across all tf.keras cell type models.
    """

    def predict(self, x, **kwarg):
        return self.training_model.predict(x)
