try:
    import tensorflow as tf
except ImportError:
    tf = None
import numpy as np


class WarmUpTraining(tf.keras.callbacks.Callback):
    def __init__(
            self,
            fn=None,
            max_beta=1.0,
            verbose=1
    ):
        super(WarmUpTraining, self).__init__()
        self.verbose = verbose
        self.beta = None
        self.fn = fn
        self.count = 0
        self.max_beta = max_beta

    def on_epoch_begin(self, epoch, logs=None):
        if epoch == 0:
            if not hasattr(self.model.loss["kl"], 'beta'):
                raise ValueError('Model must have a "beta" attribute.')

        self.beta = np.minimum(self.max_beta, (epoch + 1) / 100)   # sigmoid((x-125)/4)
        tf.keras.backend.set_value(self.model.loss["kl"].beta, self.beta)
        if self.verbose > 0:
            print('\nSet beta to %s' % self.beta)

        if self.fn is not None and epoch == 99:
            self.model.save_weights(self.fn + "_weights_wup.h5")
            if self.verbose > 0:
                print('\nSave weights after warm-up training at epoch %03d' % (epoch + 1))


class LearningRateSchedule(tf.keras.callbacks.Callback):
    def __init__(
            self,
            verbose=1
    ):
        super(LearningRateSchedule, self).__init__()
        self.verbose = verbose

    def on_epoch_end(self, epoch, logs=None):
        # pseudo_inputs = self.model.layers[3].pseudo_inputs_layer.pseudo_inputs.numpy()
        # print("\nepoch %s pseudo inputs: %s / %s / %s  (max / mean / min) " %
        #       (epoch + 1, pseudo_inputs.max(), pseudo_inputs.mean(), pseudo_inputs.min()))
        if epoch == 199:
            lr = tf.keras.backend.get_value(self.model.optimizer.lr)
            tf.keras.backend.set_value(self.model.optimizer.lr, lr / 10)
            if self.verbose > 0:
                print('\nReduce lr training at epoch %03d to %s' % (epoch + 1, lr / 10))

        if epoch == 249:
            lr = tf.keras.backend.get_value(self.model.optimizer.lr)
            tf.keras.backend.set_value(self.model.optimizer.lr, lr / 10)
            if self.verbose > 0:
                print('\nReduce lr training at epoch %03d to %s' % (epoch + 1, lr / 10))

        if epoch == 299:
            self.model.stop_training = True
            if self.verbose > 0:
                print('\nTerminate training at epoch %03d' % (epoch + 1))
