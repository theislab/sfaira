import numpy as np
try:
    import tensorflow as tf
except ImportError:
    tf = None


class LossLoglikelihoodNb(tf.keras.losses.Loss):

    def __init__(self, average=True, *args, **kwargs):
        super(LossLoglikelihoodNb, self).__init__(*args, **kwargs, name="ll_nb")
        self.average = average

    def call(
            self,
            y_true,
            y_pred
    ):
        """Implements the negative log likelihood loss as VAE reconstruction loss"""
        x = y_true
        loc, scale = tf.split(y_pred, num_or_size_splits=2, axis=1)

        eta_loc = tf.math.log(loc)
        eta_scale = tf.math.log(scale)

        log_r_plus_mu = tf.math.log(scale + loc)

        ll = tf.math.lgamma(scale + x)
        ll = ll - tf.math.lgamma(x + tf.ones_like(x))
        ll = ll - tf.math.lgamma(scale)
        ll = ll + tf.multiply(x, eta_loc - log_r_plus_mu) + tf.multiply(scale, eta_scale - log_r_plus_mu)

        ll = tf.clip_by_value(ll, -300, 300, "log_probs")
        neg_ll = -ll
        if self.average:
            neg_ll = tf.reduce_mean(neg_ll)
        else:
            # sum over features, average over batch
            neg_ll = tf.reduce_mean(tf.reduce_sum(neg_ll, axis=1), axis=0)
        return neg_ll


class LossLoglikelihoodGaussian(tf.keras.losses.Loss):

    def __init__(self, average=True, *args, **kwargs):
        super(LossLoglikelihoodGaussian, self).__init__(*args, **kwargs, name="ll_norm")
        self.average = average

    def call(
            self,
            y_true,
            y_pred
    ):
        """Implements the gaussian log likelihood loss as VAE reconstruction loss"""
        loc, scale = tf.split(y_pred, num_or_size_splits=2, axis=1)

        ll = -tf.math.log(scale * tf.math.sqrt(2. * np.pi)) - 0.5 * tf.math.square((y_true - loc) / scale)
        ll = tf.clip_by_value(ll, -300, 300, "log_probs")
        neg_ll = -ll
        if self.average:
            neg_ll = tf.reduce_mean(neg_ll)
        else:
            # sum over features, average over batch
            neg_ll = tf.reduce_mean(tf.reduce_sum(neg_ll, axis=1), axis=0)
        return neg_ll


class LossCrossentropyAgg(tf.keras.losses.Loss):

    def __init__(self, *args, **kwargs):
        super(LossCrossentropyAgg, self).__init__(*args, **kwargs, name="cce_agg")

    def call(
            self,
            y_true,
            y_pred
    ):
        """ Modified crossentropy that aggregates allowed output classes into single class. """
        y_pred = tf.clip_by_value(y_pred, 1e-10, 1., "y_pred")
        ll_cce_agg = -tf.math.log(tf.reduce_mean(y_true * y_pred, axis=1, keepdims=False))
        return ll_cce_agg


class KLLoss(tf.keras.losses.Loss):

    def __init__(self):
        super(KLLoss, self).__init__()
        self.beta = tf.Variable(1.0, dtype=tf.float32, trainable=False)

    def call(
            self,
            y_true,
            y_pred
    ):
        expected_logqz_x, expected_logpz = tf.split(y_pred, num_or_size_splits=2, axis=1)

        kl_loss = tf.reduce_mean(expected_logqz_x - expected_logpz, axis=0)
        return self.beta * kl_loss
