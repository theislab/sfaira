try:
    import tensorflow as tf
except ImportError:
    tf = None


class NegBinOutput(tf.keras.layers.Layer):
    """Negative binomial output layer"""

    def __init__(
            self,
            original_dim=None,
            name='neg_bin_output',
            **kwargs
    ):

        super().__init__(name=name, **kwargs)

        self.means = tf.keras.layers.Dense(original_dim, activation='linear')
        self.var = tf.keras.layers.Dense(original_dim, activation='linear')

    def call(self, inputs, **kwargs):
        activation, sf = inputs
        mean, var = self.means(activation), self.var(activation)

        # clip to log of largest values supported by log operation
        bound = 60.
        mean_clip = tf.clip_by_value(mean, -bound, bound, "decoder_clip")
        var_clip = tf.clip_by_value(var, -bound, bound, "decoder_clip")

        invlinker_mean = tf.exp(mean_clip + sf)
        invlinker_var = tf.exp(var_clip)

        return [invlinker_mean, invlinker_var]


class NegBinSharedDispOutput(tf.keras.layers.Layer):
    """Negative binomial output layer with a single dispersion estimate per features"""

    def __init__(
            self,
            original_dim=None,
            name='neg_bin_shared_disp_output',
            **kwargs
    ):

        super().__init__(name=name, **kwargs)

        self.means = tf.keras.layers.Dense(original_dim, activation='linear')
        self.var = self.add_weight(
            "var_bias",
            shape=[1, original_dim]
        )

    def call(self, inputs, **kwargs):
        activation, sf = inputs
        mean = self.means(activation)
        var = self.var
        var = tf.broadcast_to(var, tf.shape(mean))

        # clip to log of largest values supported by log operation
        bound = 60.
        mean_clip = tf.clip_by_value(mean, -bound, bound, "decoder_clip")
        var_clip = tf.clip_by_value(var, -bound, bound, "decoder_clip")

        invlinker_mean = tf.exp(mean_clip + sf)
        invlinker_var = tf.exp(var_clip)

        return [invlinker_mean, invlinker_var]


class NegBinConstDispOutput(tf.keras.layers.Layer):
    """Negative binomial output layer with dispersion set as constant (=1)."""

    def __init__(
            self,
            original_dim=None,
            name='neg_bin_const_disp_output',
            **kwargs
    ):

        super().__init__(name=name, **kwargs)

        self.means = tf.keras.layers.Dense(original_dim, activation='linear')
        self.var_constant = 1.

    def call(self, inputs, **kwargs):
        activation, sf = inputs
        mean = self.means(activation)
        var = tf.constant([[self.var_constant]], dtype=activation.dtype)
        var = tf.broadcast_to(var, tf.shape(mean))

        # clip to log of largest values supported by log operation
        bound = 60.
        mean_clip = tf.clip_by_value(mean, -bound, bound, "decoder_clip")
        var_clip = tf.clip_by_value(var, -bound, bound, "decoder_clip")

        invlinker_mean = tf.exp(mean_clip + sf)
        invlinker_var = tf.exp(var_clip)

        return [invlinker_mean, invlinker_var]


class GaussianOutput(tf.keras.layers.Layer):
    """
    Gaussian output layer.

    Size factor only makes sense if logged and data is positive and logged.
    """

    def __init__(
            self,
            original_dim=None,
            name='gaussian_output',
            **kwargs
    ):

        super().__init__(name=name, **kwargs)

        self.means = tf.keras.layers.Dense(original_dim, activation='linear')
        self.var = tf.keras.layers.Dense(original_dim, activation='linear')

    def call(self, inputs, **kwargs):
        activation, sf = inputs
        mean, var = self.means(activation), self.var(activation)

        # clip to log of largest values supported by log operation
        bound = 60.
        mean_clip = tf.clip_by_value(mean, tf.exp(-bound), tf.exp(bound), "decoder_clip")
        var_clip = tf.clip_by_value(var, -bound, bound, "decoder_clip")

        invlinker_mean = mean_clip + sf
        invlinker_var = tf.exp(var_clip)

        return [invlinker_mean, invlinker_var]


class GaussianSharedStdOutput(tf.keras.layers.Layer):
    """
    Gaussian output layer with a single standard deviation estimate per features.

    Size factor only makes sense if logged and data is positive and logged.
    """

    def __init__(
            self,
            original_dim=None,
            name='gaussian_shared_disp_output',
            **kwargs
    ):

        super().__init__(name=name, **kwargs)

        self.means = tf.keras.layers.Dense(original_dim, activation='linear')
        self.var = self.add_weight(
            "var_bias",
            shape=[1, original_dim]
        )

    def call(self, inputs, **kwargs):
        activation, sf = inputs
        mean = self.means(activation)
        var = self.var
        var = tf.broadcast_to(var, tf.shape(mean))

        # clip to log of largest values supported by log operation
        bound = 60.
        mean_clip = tf.clip_by_value(mean, tf.exp(-bound), tf.exp(bound), "decoder_clip")
        var_clip = tf.clip_by_value(var, -bound, bound, "decoder_clip")

        invlinker_mean = mean_clip + sf
        invlinker_var = tf.exp(var_clip)

        return [invlinker_mean, invlinker_var]


class GaussianConstStdOutput(tf.keras.layers.Layer):
    """
    Gaussian output layer with standard deviation set as constant (=1).

    Size factor only makes sense if logged and data is positive and logged.
    """

    def __init__(
            self,
            original_dim=None,
            name='gaussian_const_disp_output',
            **kwargs
    ):

        super().__init__(name=name, **kwargs)

        self.means = tf.keras.layers.Dense(original_dim, activation='linear')
        self.var_constant = 1.

    def call(self, inputs, **kwargs):
        activation, sf = inputs
        mean = self.means(activation)
        var = tf.constant([[self.var_constant]], dtype=activation.dtype)
        var = tf.broadcast_to(var, tf.shape(mean))

        # clip to log of largest values supported by log operation
        bound = 60.
        mean_clip = tf.clip_by_value(mean, tf.exp(-bound), tf.exp(bound), "decoder_clip")
        var_clip = tf.clip_by_value(var, -bound, bound, "decoder_clip")

        invlinker_mean = mean_clip + sf
        invlinker_var = tf.exp(var_clip)

        return [invlinker_mean, invlinker_var]
