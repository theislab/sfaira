try:
    import tensorflow as tf
except ImportError:
    tf = None


class PreprocInput(tf.keras.layers.Layer):

    def __init__(
            self,
            name='preproc_input',
            **kwargs
    ):
        super().__init__(name=name, **kwargs)

    def call(self, inputs, **kwargs):
        sf = tf.reduce_sum(inputs, axis=1, keepdims=True) / 10000.
        inputs = tf.divide(inputs, sf)
        inputs = tf.math.log(inputs + 1)
        return inputs
