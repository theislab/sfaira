import numpy as np
try:
    import tensorflow as tf
except ImportError:
    tf = None


def custom_mse(y_true, y_pred, sample_weight=None):
    y_pred = tf.split(y_pred, num_or_size_splits=2, axis=1)[0]
    se = tf.square(tf.subtract(y_true, y_pred))
    se_red = tf.reduce_mean(se)
    return se_red


def custom_mean_squared_logp1_error(y_true, y_pred, sample_weight=None):
    y_pred = tf.split(y_pred, num_or_size_splits=2, axis=1)[0]
    y_true = tf.math.log(y_true + 1.)
    y_pred = tf.math.log(y_pred + 1.)
    se = tf.square(tf.subtract(y_true, y_pred))
    se_red = tf.reduce_mean(se)
    return se_red


def custom_negll_nb(y_true, y_pred, sample_weight=None):
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
    neg_ll = tf.reduce_mean(neg_ll)
    return neg_ll


def custom_negll_gaussian(y_true, y_pred, sample_weight=None):
    loc, scale = tf.split(y_pred, num_or_size_splits=2, axis=1)

    ll = -tf.math.log(scale * tf.math.sqrt(2. * np.pi)) - 0.5 * tf.math.square((y_true - loc) / scale)
    ll = tf.clip_by_value(ll, -300, 300, "log_probs")
    neg_ll = -ll
    neg_ll = tf.reduce_mean(neg_ll)
    return neg_ll


def custom_kl(y_true, y_pred, sample_weight=None):
    expected_logqz_x, expected_logpz = tf.split(y_pred, num_or_size_splits=2, axis=1)
    kl_loss = tf.reduce_mean(expected_logqz_x - expected_logpz)
    return kl_loss


def custom_cce_agg(y_true, y_pred, sample_weight=None):
    y_pred = tf.clip_by_value(y_pred, 1e-10, 1., "y_pred")
    ll_cce_agg = -tf.math.log(tf.reduce_mean(y_true * y_pred, axis=1, keepdims=False))
    return ll_cce_agg


class CustomAccAgg(tf.keras.metrics.Metric):

    def __init__(self, name='acc_agg', **kwargs):
        super(CustomAccAgg, self).__init__(name=name, **kwargs)
        self.acc_agg = self.add_weight(name='acc', initializer='zeros')
        self.count = self.add_weight(name='count', initializer='zeros')

    def update_state(self, y_true, y_pred, sample_weight=None):
        phat_pos_agg = tf.reduce_sum(y_true * y_pred, axis=1, keepdims=True)
        acc_agg = tf.cast(
            phat_pos_agg > tf.reduce_max((tf.ones_like(y_true) - y_true) * y_pred, axis=1),
            dtype=y_true.dtype
        )
        # Do not use weighting for accuracy.
        self.acc_agg.assign_add(tf.reduce_mean(acc_agg))
        self.count.assign_add(1.)

    def result(self):
        return tf.divide(self.acc_agg, self.count)

    def reset_states(self):
        self.acc_agg.assign(0.)
        self.count.assign(0.)


class CustomTprClasswise(tf.keras.metrics.Metric):

    def __init__(self, k: int, name='tpr', **kwargs):
        super(CustomTprClasswise, self).__init__(name=name, **kwargs)
        self.tp = self.add_weight(shape=(k,), name='tp', initializer='zeros')
        self.fn = self.add_weight(shape=(k,), name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred, sample_weight=None):
        tp_by_class = tf.reduce_sum(tf.cast(
            y_pred == tf.reduce_max(y_pred, axis=1, keepdims=True),
            dtype=y_true.dtype
        ) * y_true, axis=0)
        fn_by_class = tf.reduce_sum(tf.cast(
            y_pred < tf.reduce_max(y_pred, axis=1, keepdims=True),
            dtype=y_true.dtype
        ) * y_true, axis=0)
        self.tp.assign_add(tp_by_class)
        self.fn.assign_add(fn_by_class)

    def result(self):
        # Catch division by zero, in that case tpr becomes zero after clipping the divisor to 1.
        divisor = tf.clip_by_value(self.tp + self.fn, 1, np.inf, "divisor")
        tpr = tf.reduce_mean(self.tp / divisor)
        return tpr

    def reset_states(self):
        self.tp.assign(tf.zeros_like(self.tp))
        self.fn.assign(tf.zeros_like(self.fn))


class CustomFprClasswise(tf.keras.metrics.Metric):

    def __init__(self, k: int, name='fpr', **kwargs):
        super(CustomFprClasswise, self).__init__(name=name, **kwargs)
        self.fp = self.add_weight(shape=(k,), name='fp', initializer='zeros')
        self.tn = self.add_weight(shape=(k,), name='tn', initializer='zeros')

    def update_state(self, y_true, y_pred, sample_weight=None):
        fp_by_class = tf.reduce_sum(tf.cast(
            y_pred == tf.reduce_max(y_pred, axis=1, keepdims=True),
            dtype=y_true.dtype
        ) * (1. - y_true), axis=0)
        tn_by_class = tf.reduce_sum(tf.cast(
            y_pred < tf.reduce_max(y_pred, axis=1, keepdims=True),
            dtype=y_true.dtype
        ) * (1. - y_true), axis=0)
        self.fp.assign_add(fp_by_class)
        self.tn.assign_add(tn_by_class)

    def result(self):
        # Catch division by zero, in that case fpr becomes zero after clipping the divisor to 1.
        divisor = tf.clip_by_value(self.fp + self.tn, 1, np.inf, "divisor")
        fpr = tf.reduce_mean(self.fp / divisor)
        return fpr

    def reset_states(self):
        self.fp.assign(tf.zeros_like(self.fp))
        self.tn.assign(tf.zeros_like(self.tn))


class CustomF1Classwise(tf.keras.metrics.Metric):

    def __init__(self, k: int, name='f1', **kwargs):
        super(CustomF1Classwise, self).__init__(name=name, **kwargs)
        self.tp = self.add_weight(shape=(k,), name='tp', initializer='zeros')
        self.fp = self.add_weight(shape=(k,), name='fp', initializer='zeros')
        self.fn = self.add_weight(shape=(k,), name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred, sample_weight=None):
        tp_by_class = tf.reduce_sum(tf.cast(
            y_pred == tf.reduce_max(y_pred, axis=1, keepdims=True),
            dtype=y_true.dtype
        ) * y_true, axis=0)
        fp_by_class = tf.reduce_sum(tf.cast(
            y_pred == tf.reduce_max(y_pred, axis=1, keepdims=True),
            dtype=y_true.dtype
        ) * (1. - y_true), axis=0)
        fn_by_class = tf.reduce_sum(tf.cast(
            y_pred < tf.reduce_max(y_pred, axis=1, keepdims=True),
            dtype=y_true.dtype
        ) * y_true, axis=0)
        self.tp.assign_add(tp_by_class)
        self.fp.assign_add(fp_by_class)
        self.fn.assign_add(fn_by_class)

    def result(self):
        # Catch divisions by zero, in that case precision or recall become zero after clipping the divisor to 1.
        divisor_precision = tf.clip_by_value(self.tp + self.fp, 1, np.inf, "divisor")
        divisor_recall = tf.clip_by_value(self.tp + self.fn, 1, np.inf, "divisor")
        precision = self.tp / divisor_precision
        recall = self.tp / divisor_recall
        precision = tf.clip_by_value(precision, 1e-100, np.inf, "divisor")
        recall = tf.clip_by_value(recall, 1e-100, np.inf, "divisor")
        f1 = tf.reduce_mean(2 * 1 / (1 / precision + 1 / recall))
        return f1

    def reset_states(self):
        self.tp.assign(tf.zeros_like(self.tp))
        self.fp.assign(tf.zeros_like(self.fp))
        self.fn.assign(tf.zeros_like(self.fn))
