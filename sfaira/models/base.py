import abc
try:
    import tensorflow as tf
except ImportError:
    tf = None


class BasicModel(abc.ABC):
    """
    This base class defines model attributes shared across all models.
    """
    _version: str
    _topology_id: str
    genome_size: int
    hyperparam: dict
    model_class: str
    model_type: str

    @property
    def version(self):
        return self._version


class BasicModelKeras(BasicModel):
    """
    This base class defines model attributes shared across all tf.keras models.
    """

    training_model: tf.keras.Model
